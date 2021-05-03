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
	"os"
	"runtime"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
)

//Global variables... Sometimes, you gotta use'em
var verb int

//If v is true, prints the d arguments to stderr
//otherwise, does nothing.
func LogV(v int, vref int, d ...interface{}) {
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
	tempinterval := flag.Float64("tempInterval", 5.0, "the interval of temperature between replicas.")
	multi := flag.Int("multi", 1, "multiplicity of the system")
	charge := flag.Int("charge", 0, "charge for the system")
	verbose := flag.Int("verbose", 0, "Level of verbosity, the higher, the more verbose.")
	exrate := flag.Int("exchangerate", 2, "Every how many ps should replica exchanges be attempted?")
	Tc := flag.Float64("tlow", 310.15, "The temperature of the coldest replica")
	Th := flag.Float64("thot", 350.15, "The temperature of the hotest replica")
	dielectric := flag.Float64("dielectric", 80.0, "The dielectric constant for continuum solvent in QM calculations. Only some values are allowed (see code) -1 for vacuum calculations. Default is acetone")
	flag.Parse()
	verb = *verbose
	//just in case.
	args := flag.Args()
	geoname := args[0]
	totaltime, err := strconv.Atoi(args[1])
	CErr(err, "main")
	fmt.Printf("Use:\n  $REPATH/xtbRE [FLAGS] geometry mdtime \n")
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
	temps := make([]float64, 0, (int(*Th)-int(*Tc))/int(*tempinterval)+2)
	for i := 0.0; i+*Tc <= *Th; i += *tempinterval {
		temps = append(temps, *Tc+i)
	}
	if *cpus < 0 {
		*cpus = runtime.NumCPU()
	}

	MDs(mol, *exrate, totaltime, *method, temps, *dielectric, *cpus)
}
