//
//    MCL: MiMiC Communication Library
//    Copyright (C) 2015-2017  Viacheslav Bolnykh,
//                             Jógvan Magnus Haugaard Olsen,
//                             Simone Meloni,
//                             Emiliano Ippoliti,
//                             Paolo Carloni
//                             and Ursula Röthlisberger.
//
//    This file is part of MCL.
//
//    MCL is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as
//    published by the Free Software Foundation, either version 3 of
//    the License, or (at your option) any later version.
//
//    MCL is distributed in the hope that it will be useful, but
//    WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <string>
#include "BaseMessageData.h"

#ifndef MIMICCOMMLIB_SYSTEMDATA_H
#define MIMICCOMMLIB_SYSTEMDATA_H

/**
 * Class encapsulatin the data about the system of the client code
 */
struct SystemData : public BaseMessageData {

public:
    SystemData() : BaseMessageData() { }

    SystemData(int type, int entity_number) : BaseMessageData(type, entity_number) {
    }

    /**
     * Atom indices
     */
    std::vector<int> indices;

    /**
     * Atom types
     */
    std::vector<std::string> types;

    /**
     * Coordinates of items mapped onto 1D array
     */
    std::vector<double> coordinates;

    /**
     * Atom masses
     */
    std::vector<double> masses;

    /**
     * Maximum order of multipoles that is used by the code
     */
    int multipole_order;

};

#endif //MIMICCOMMLIB_SYSTEMDATA_H
