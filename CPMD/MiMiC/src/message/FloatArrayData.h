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

#include "BaseMessageData.h"
#include <vector>

#ifndef MIMICCOMMLIB_FLOATARRAYDATA_H
#define MIMICCOMMLIB_FLOATARRAYDATA_H

/**
 * Class for sending an array of doubles (like coordinates)
 */
struct FloatArrayData : public BaseMessageData {
public:
    /**
     * Array of double values to be sent
     */
    std::vector<double> array;

    FloatArrayData() : BaseMessageData() { }

    FloatArrayData(int type, int entity_number) : BaseMessageData(type, entity_number) {
    }

};
#endif //MIMICCOMMLIB_FLOATARRAYDATA_H
