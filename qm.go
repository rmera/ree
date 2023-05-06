/*
 * qm.go, part of xtbRE
 *
 *
 * Copyright 2021 Raul Mera <rmera{at}usach(dot)cl>
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
 */

package main

import (
	"bufio"
	"fmt"
	"math"
	rand "math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/qm"
)

//This has been tested and seems to work fine.

type M struct {
	T float64
	V float64
}

// A small function to log the temperatures by which our "Reference" system goes.
//It will print a message at verbose level "2" if a switch happens involving the reference temperature.
func SignalT(T1, T2, Tref float64) {
	if T1 == Tref {
		LogV(2, "T exchanged:", T2)
	} else if T2 == Tref {
		LogV(2, "T exchanged:", T1)
	}
}

//This is the "central control" that holds the info for all the replicas, runs them
//collect their outputs, decides whether each pair should switch temperatures, and, if so,
//sends each its new T. It also deals with the action needed when all replicas finish running.
func MDs(mol *chem.Molecule, ctime, ttime int, method string, temps []float64, dielectric float64, cpus int) {
	//hopefully this seed is good enough.
	rand.Seed(time.Now().UTC().UnixNano())
	Tref := temps[0] //This is the reference T, the simulation at this T is the one we actually care about.
	rcpus := cpus / len(temps)
	if rcpus < 1 {
		rcpus = 1
	}
	Hslice := make([]*RHandler, 0, len(temps))
	for i, v := range temps {
		h := NewRHandler(i, v)
		Hslice = append(Hslice, h)
	}
	//Now we prepare an launch the replicas. Scary stuff
	Hs := Handlers(Hslice)
	for _, v := range Hs {
		R := &Replica{ID: v.ID, T: v.CurrT, Time: ctime, TotTime: ttime, epsilon: dielectric, C: v.C, Tref: Tref}
		go RunReplica(R, method, mol, rcpus)
	}
	//OK, now, way scarier, we need to goo on collecting data for all molecules in a loop, until we get "finishing" values
	//the preopt
	Q := new(qm.Calc)
	Q.Method = method
	Q.Job = qm.Job{Opti: true}
	Q.Dielectric = dielectric
	xtb := qm.NewXTBHandle()
	xtb.SetnCPU(cpus)
	xtb.BuildInput(mol.Coords[0], mol, Q)
	xtb.Run(true)
	var err error
	mol.Coords[0], err = xtb.OptimizedGeometry(mol)
	CErr(err, "Initial optimization failed")
	NewEven := func(even bool) bool {
		return !even
	}
	var even bool = true
	var term bool
	for {
		//Note that we don't use "select"; we wait until all workers are
		//ready before starting the next cycle. This is how we ensure synchronization.
		for _, v := range Hs {
			tv := <-v.C.repdata
			v.CurrT = tv[0]
			v.CurrP = tv[1]
			//a Temperature <=0 signals that
			//workers reached the end of the total simulation time.
			if v.CurrT <= 0 {
				term = true
			}

		}

		if term {
			for _, v := range Hs {
				v.C.Close()
			}
			break
		}
		//From now on, Hs are sorted by temperature
		sort.Sort(Hs)
		//we iterate the list from higher to lowe temperature
		//we need i>0 because we will be
		//using both i and i-1

		for i := len(Hs) - 1; i > 0; i-- {

			//In the new system, even i's exchange on even cycles
			//odd i's exchange on odd cycles. NewEven is the function
			//that ensures that "even" is true on even cycles and false
			//on odd cycles. With this system, we try to approach
			//Gromac's behaviour.
			if i%2 != 0 && even {
				continue
			} else if i%2 == 0 && !even {
				continue
			}

			/**********This is the old system **********/
			//when we get to the "i" replica (save for i==len(Hs)-1) the replica could have already
			//exchanged T with the i+1 one. If that happened (which Hs[i].Switched marks), we do not
			//attempt to exchange it with i-1.
			//if Hs[i].Switched || Hs[i].SwitchedPrev {
			//	continue
			//}
			/*******End old system********************/
			p := Metropolis(&M{T: Hs[i].CurrT, V: Hs[i].CurrP}, &M{T: Hs[i-1].CurrT, V: Hs[i-1].CurrP})
			if p == 1 || rand.Float64() < p {
				Hs[i].C.newT <- Hs[i-1].CurrT
				Hs[i-1].C.newT <- Hs[i].CurrT
				Hs[i-1].Switched = true //these things are remanants of the old system. I'll just leave them for now.
				Hs[i].SwitchedPrev = true
				SignalT(Hs[i].CurrT, Hs[i-1].CurrT, Tref)

			} else {
				Hs[i].C.newT <- Hs[i].CurrT

			}
			even = NewEven(even) //if even was true, it becomes false, and vice versa
		}
		//since the previous look only reaches down to Hs[1], if Hs[1] doesn't switch with Hs[0], Hs[0] will not get the "continue" signal
		//and will get stuck, so here we send the signal by hand if that happened.
		if !Hs[0].Switched {
			Hs[0].C.newT <- Hs[0].CurrT
		}
		//We need to reset the "Switched" now, so we can start the next cycle
		for _, v := range Hs {
			v.Switched = false
		}
		//	panic("V") ///////

	}
	//Now we need to stitch together the different parts of the trajectory at the lowest T, which will be our final traj.
	cq := exec.Command("sh", "-c", "cat *.trj > fulltraj.xyz")
	cq.Run()

}

//Implements the Metropolis function (or so I hope)
//This is where I have the most doubts.
func Metropolis(i, j *M) float64 {
	b1 := 1 / (chem.R * i.T)
	b2 := 1 / (chem.R * j.T)
	Delta := (b1 - b2) * (j.V - i.V)
	if Delta <= 0 {
		return 1.0
	} else {
		return math.Exp(-1 * Delta)
	}

}

//This structure contain the two channels each replica has to communicate with the "control center"
//newT where it gets the new temperature from "control" and repdata, where it sends it data to control.
type com struct {
	newT    chan float64   //The program will send the new temperature to the gorutine
	repdata chan []float64 // Signals the gorutine is waiting for information to arrive on the T channel. A "false" value means the gorutine has finished (and the channel, closed).
}

func (C *com) Close() {
	close(C.newT)
	close(C.repdata)

}

//when a replica is ready, it will send its last temperature
//to the 'center' via oldT.
//the 'center' will use the last temperatures from workers to
//decide which replicas are exchanging temperatures, and then
//send the new temperatures to each via newT
func Newcom() *com {
	c := new(com)
	//they need to be buffered
	c.newT = make(chan float64, 1)
	c.repdata = make(chan []float64, 1)
	return c
}

//This structure  has the "control side" info for each replica.
//The C set of channels is shared with the "replica side"
//structure, so control and replicas can communicate.
type RHandler struct {
	ID           int
	C            *com
	CurrT        float64
	CurrP        float64
	Switched     bool
	SwitchedPrev bool
}

func NewRHandler(ID int, initemp float64) *RHandler {
	H := new(RHandler)
	H.C = Newcom()
	H.CurrT = initemp
	H.ID = ID
	return H
}

type Handlers []*RHandler

func (H Handlers) Less(i, j int) bool {
	return H[i].CurrT < H[i].CurrT
}

func (H Handlers) Len() int {
	return len(H)
}

func (H Handlers) Swap(i, j int) {
	H[i], H[j] = H[j], H[i]
}

//This is the "replica side" information for each replica
//the C set of channels is shared with the corresponding "control side" RHandler
//structure, so both sides can communicate.
type Replica struct {
	ID      int
	T       float64 //The current temperature
	Tref    float64 //the reference temperature for all replicas (the lowest temperature)
	PrevTs  []float64
	V       []float64 //potential energies
	TotTime int       //ps
	epsilon float64
	Time    int //ps
	C       *com
}

//Run replica controls the run of one of the replicas in the system. It takes a *Replica structure with the data needed for the whole thing
//plus the geometry (for reading only!) and the CPUs to be used in the calculation.
//it will run simulations of R.Time lenght until a total R.TotTime simulation time is reached. After each simulation run it will signal to
//the "center" the temperature it has been using, and wait for a new temperature from the center. It will set  the temperature to the received value
//scale the velocities by sqrt(newT/oldT) and restart the simulation. When it has finished the total sim. time. it will report -1 instead of a temperature,
//and close the channels. RunReplica will also keep track of the number of simulations ran, and all the temperatures used. It will save all trajectorie
//identified with their number and temperature (so the program can later build a whole trajectory at a given temperature with all the pieces from different
//replicas, by simply concatenating them together.
func RunReplica(R *Replica, method string, ReadMol *chem.Molecule, cpus int) {
	dirname := strconv.Itoa(R.ID)
	err := os.Mkdir(dirname, os.FileMode(0755))                      //we don't need to catch this error, we'll just get int in the next statement.
	if err != nil && !strings.Contains(err.Error(), "file exists") { //If the dir was already there, we just continue.
		CErr(err, "worker "+strconv.Itoa(R.ID))
	}
	Q := new(qm.Calc)
	Q.Method = method
	Q.Job = qm.Job{MD: true}
	Q.MDTime = R.Time
	Q.MDTemp = R.T
	Q.Dielectric = R.epsilon
	xtb := qm.NewXTBHandle()
	xtb.SetnCPU(cpus)
	xtb.SetWorkDir(dirname)
	xtb.BuildInput(ReadMol.Coords[0], ReadMol, Q)
	//the input for the single point1
	Qsp := new(qm.Calc)
	Qsp.Method = method
	Qsp.Job = qm.Job{SP: true}
	xtbsp := qm.NewXTBHandle()
	xtbsp.SetnCPU(cpus)
	xtbsp.SetName("SP")
	xtbsp.SetWorkDir(dirname)
	xtbsp.BuildInput(ReadMol.Coords[0], ReadMol, Qsp)
	elapsed := 0
	R.PrevTs = append(R.PrevTs, R.T)
	for {
		xtb.Run(true)
		lastT := R.PrevTs[len(R.PrevTs)-1]
		elapsed += R.Time
		if lastT == R.Tref {
			os.Rename(dirname+"/xtb.trj", fmt.Sprintf("xtb_%08d_%5.2f_.trj", len(R.PrevTs), lastT))
			weaita, _ := os.Create(fmt.Sprintf("%s/%06d", dirname, len(R.PrevTs)))
			weaita.Close()
			//	cq := exec.Command("sh", "-c", fmt.Sprintf("cp xtb.trj xtb_%d_%5.2f_.trj", len(R.PrevTs), lastT)) //yeah, I'm lazy. Also this is Unix-only program anyway

			//	cq.Dir = dirname
			//	cq.Run()
		} else {
			os.Remove(dirname + "/xtb.trj")
		}
		//We will try to remove scoord files left by xtb, but if it doesn't work, it doesn't work
		//the program will just keep running.
		toremove, err := filepath.Glob(dirname + "/scoord*")
		if err == nil {
			for _, f := range toremove {
				_ = os.Remove(f)
			}
		}

		//	cq := exec.Command("sh", "-c", fmt.Sprintf("cp mdrestart mdrestart_%d_%5.2f", len(R.PrevTs), lastT)) ///////
		//	cq.Dir = dirname                                                                                     ////////////
		//	cq.Run()                                                                                             //////////
		if elapsed >= R.TotTime {
			R.C.repdata <- []float64{-1, -1}
			//R.C.Close()
			break
		}
		pot, err := lastPot(xtbsp, ReadMol, dirname)
		CErr(err, fmt.Sprintf("Worker %d, while obtaining the last potential energy of chunk %d from %5.3f", R.ID, len(R.PrevTs), lastT))
		R.C.repdata <- []float64{lastT, pot}
		newT := <-R.C.newT
		//Maybe this is too dirty. When the program actually works I'll consider switching this for a proper Go function.
		//lol who am I kidding.
		cq := exec.Command("sh", "-c", fmt.Sprintf("sed -i s/temp=%5.2f/temp=%5.2f/g %s/gochem.inp", lastT, newT, dirname))
		cq.Run()
		cq = exec.Command("sh", "-c", fmt.Sprintf("sed -i s/restart=false/restart=true/g %s/gochem.inp", dirname)) //if restart is already true this will do nothing
		cq.Run()
		err = ScaleVel(lastT, newT, ReadMol.Len(), dirname)
		CErr(err, fmt.Sprintf("Worker %d, while scaling chunk %d from %5.3f to %5.3f", R.ID, len(R.PrevTs), lastT, newT))
		//		cq = exec.Command("sh", "-c", fmt.Sprintf("cp mdrestart mdrestartRM_%d_%5.2f", len(R.PrevTs), lastT)) ///////
		//		cq.Dir = dirname                                                                                      ////////////
		//		cq.Run()                                                                                              //////////

		R.PrevTs = append(R.PrevTs, newT)

	}

}

//I honestly don't quite know how to ensure that I get the potential energy of the last geometry,
//other than this, to just read the geometry from the mdrestart file, and run an SP calculation.
func lastPot(xtbsp *qm.XTBHandle, molR *chem.Molecule, dirname string) (float64, error) {
	err := lastGeo(molR, dirname)
	if err != nil {
		return 0, err
	}
	xtbsp.Run(true)
	return xtbsp.Energy()

}

//Creates a "SP.xyz"  file with the last coordinates from an xtb MD simulation
//in the current directory
func lastGeo(mol *chem.Molecule, dirname string) error {
	pos := make([][]float64, 0, mol.Len())
	fin, err := os.Open(dirname + "/mdrestart")
	defer fin.Close() //If no errors, we'll close fin manually before the function exists. That's fine, we won't check for any errors.
	if err != nil {
		return err
	}
	bfin := bufio.NewReader(fin)
	_, err = bfin.ReadString('\n') //the first line doesn't have coordinates
	if err != nil {
		return err
	}
	for i := 0; i < mol.Len(); i++ {
		line, err := bfin.ReadString('\n')
		if err != nil {
			return err //I *think* I should never find an EOF since I know the file should be mol.Len() lines, so this error should be always bad.
		}
		p, err := scaleline(line, 1)
		if err != nil {
			return err
		}
		pos = append(pos, p)

	}
	//I assume the atomic positions are in Bohrs
	lg, err := os.Create(dirname + "/SP.xyz")
	if err != nil {
		return err
	}
	b2a := chem.Bohr2A
	lg.WriteString(fmt.Sprintf("%d\n\n", mol.Len()))
	for i, v := range pos {
		at := mol.Atom(i)
		str := fmt.Sprintf("%s  %-8.5f  %-8.5f  %-8.5f\n", at.Symbol, v[0]*b2a, v[1]*b2a, v[2]*b2a)
		_, err = lg.WriteString(str)
		if err != nil {
			return err
		}
	}
	return nil

}

func scaleline(s string, scalefac float64) ([]float64, error) {
	sf := make([]string, 3, 3)
	sf[0] = strings.Replace(s[0:22], "D", "E", 1)
	sf[1] = strings.Replace(s[22:44], "D", "E", 1)
	sf[2] = strings.Replace(s[44:66], "D", "E", 1)
	ret := make([]float64, 3, 3)
	var err error
	for i := 0; i < 3; i++ {
		ns := strings.TrimSpace(sf[i])
		ret[i], err = strconv.ParseFloat(ns, 64)
		if err != nil {
			return nil, err
		}
		ret[i] = ret[i] * scalefac
	}
	return ret, nil

}

//scales all velocities by Tnew/Told as required by the RE protocol.
//note that it is the modulus of the velofity that needs to be scaled!
func ScaleVel(oldT, newT float64, mollen int, dirname string) error {
	scalefac := math.Sqrt(newT / oldT)
	//	fmt.Println("will scale", scalefac) //////////////
	pos := make([]string, 0, mollen)
	vels := make([][]float64, 0, mollen)
	//First we read/parse the original data from the mdrestart file.
	fin, err := os.Open(dirname + "/mdrestart")
	defer fin.Close() //If no errors, we'll close fin manually before the function exists. That's fine, we won't check for any errors.
	if err != nil {
		return err
	}
	bfin := bufio.NewReader(fin)
	firstline, err := bfin.ReadString('\n')
	if err != nil {
		return err //I *think* I should never find an EOF since I know the file should be mollen lines, so this error should be always bad.
	}

	for i := 0; i < mollen; i++ {
		line, err := bfin.ReadString('\n')
		if err != nil {
			return err //I *think* I should never find an EOF since I know the file should be mollen lines, so this error should be always bad.
		}
		pos = append(pos, line[0:66])
		v, err := scaleline(line[66:], scalefac)
		if err != nil {
			return err
		}
		vels = append(vels, v)

	}
	fin.Close()
	//now we write our own mdrestart with the new values
	fout, err := os.Create(dirname + "/mdrestart")
	fout.WriteString(firstline)
	defer fout.Close()
	if err != nil {
		return err
	}
	for i, v := range vels {
		str := fmt.Sprintf("%s%22.14E%22.14E%22.14E\n", pos[i], v[0], v[1], v[2])
		str = strings.Replace(str, "E", "D", -1) //xtb wants "D" instead of "E" in scientific notation.
		_, err = fout.WriteString(str)
		if err != nil {
			return err
		}
	}
	return nil
}
