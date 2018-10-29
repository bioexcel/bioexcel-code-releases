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

#ifndef MIMICCOMMLIB_MESSAGEDATA_H
#define MIMICCOMMLIB_MESSAGEDATA_H

/**
 * Abstract class for data contained in the message
 */
struct BaseMessageData {
public:
    /**
     * Type of the message (probably not needed anymore)
     */
    int type;

    /**
     * Number of entities contained in this chunk of data (will be needed if we are forced to work with c++98)
     */
    int entity_number;

    BaseMessageData() { }

    BaseMessageData(int type, int entity_number) : type(type), entity_number(entity_number) { }
    virtual ~BaseMessageData(){};

};
#endif //MIMICCOMMLIB_MESSAGEDATA_H
