/*
 * main.go, part of xtbRE
 *
 *
 * Copyright 2020 Raul Mera <rmera{at}usach(dot)cl>
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

/*To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche*/

package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
)

//Global variables... Sometimes, you gotta use'em
var vref int

//If v is true, prints the d arguments to stderr
//otherwise, does nothing.
func LogV(v int, d ...interface{}) {
	if v >= vref {
		fmt.Fprintln(os.Stderr, d...)
	}

}

func CErr(err error, info string) {
	if err != nil {
		log.Fatal(err, info)
	}
}

func main() {
	//There will be _tons_ of flags, but they are meant not to be needed the 99% of the time.
	method := flag.String("method", "gfnff", "the xTB method for the simulation")
	cpus := flag.Int("cpus", -1, "the total CPUs used for the QM calculations. If a number <0 is given, all logical CPUs are used")
	nreps := flag.Int("replicas", 10, "Maximum amount of replicas to be used. It might not be reached") //this default is for testing, it should probably be changed for production.
	multi := flag.Int("multi", 1, "multiplicity of the system")
	keeptrj := flag.Bool("keeptrj", false, "Keep the trjs for the reference temperature in the main directory")
	charge := flag.Int("charge", 0, "charge for the system")
	verbose := flag.Int("verbose", 1, "Level of verbosity, the higher, the more verbose.")
	exrate := flag.Int("exchangerate", 2, "Every how many ps should replica exchanges be attempted?")
	Tc := flag.Float64("tref", 310.15, "The temperature of the 'reference' or coldest replica")
	Th := flag.Float64("maxt", 453.15, "The highest temperature allowed for a replica. It might not be reached.")
	dielectric := flag.Float64("dielectric", 80.0, "The dielectric constant for continuum solvent in QM calculations. Only some values are allowed (see code) -1 for vacuum calculations. Default is acetone")
	flag.Parse()
	fmt.Println("REE: Performs a replica exchange simulation with xtb")
	fmt.Println("Use: ree OPTIONS geometry.xyz SIMULATION_TIME")
	fmt.Println("Use ree --help to see the available options")
	fmt.Println("The total replicas employed will depend on the")
	fmt.Println("-replicas and the -maxt options, which ever")
	fmt.Println("is reached first.")
	vref = *verbose
	//just in case.
	args := flag.Args()
	geoname := args[0]
	totaltime, err := strconv.Atoi(args[1])
	CErr(err, "main")
	//fmt.Printf("Use:\n  $REPATH/xtbRE [FLAGS] geometry mdtime \n")
	var mol *chem.Molecule
	extension := strings.ToLower(strings.Split(geoname, ".")[1])
	switch extension {
	case "gro":
		mol, err = chem.GroFileRead(geoname)
	case "pdb":
		mol, err = chem.PDBFileRead(geoname, false)
	default:
		mol, err = chem.XYZFileRead(geoname)
	}
	CErr(err, "main")
	mol.SetCharge(*charge) //needed for the MD and the partial charges calculation
	mol.SetMulti(*multi)
	temps := GenTemps(*Tc, *Th, *nreps, mol.Len())
	LogV(1, "Temperatures:", temps)
	//	temps := make([]float64, 0, (int(*Th)-int(*Tc))/int(*tempinterval)+2)
	//	for i := 0.0; i+*Tc <= *Th; i += *tempinterval {
	//		temps = append(temps, *Tc+i)
	//	}
	if *cpus < 0 {
		*cpus = runtime.NumCPU()
	}
	//Here is where the magic happens :-)
	MDs(mol, *exrate, totaltime, *method, temps, *dielectric, *cpus)

	//We have accumulated several "trj" files, which we already collected in one trajectory.
	//We delete them now, unless the user says not to.
	if !*keeptrj {
		toremove, err := filepath.Glob("*.trj")
		if err == nil {
			for _, f := range toremove {
				_ = os.Remove(f)
			}
		}
	}

}

//Gen temps generates nreps-1 temperatures for the RE procedure
//given the first one, Tc, and the number of atoms in the system.
//It uses the procedure outlined in the Gromacs 2021.1 Manual, p340
func GenTemps(Tc, Th float64, nreps, natoms int) []float64 {
	temps := make([]float64, 1, nreps)
	temps[0] = Tc
	epsilon := 1 / math.Sqrt(2*float64(natoms))
	for i := 1; i < nreps; i++ {
		T := (1 + epsilon) * temps[i-1]
		if T > Th {
			return temps
		}
		temps = append(temps, T)
	}
	return temps

}
