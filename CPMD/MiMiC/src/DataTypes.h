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

#ifndef MIMICCOMMLIB_DATATYPES_H
#define MIMICCOMMLIB_DATATYPES_H

/**
 * \file DataTypes.h
 * \brief Enum containing data types supported by the library
 */
enum DataType{

    /**
     * \var double type
     */
    TYPE_DOUBLE = 1,

    /**
     * \var integer type
     */
    TYPE_INT = 2,

    /**
     * \var character type
     */
    TYPE_CHAR = 3,

    /**
     * \var float type
     */
    TYPE_FLOAT = 4
};
#endif //MIMICCOMMLIB_DATATYPES_H
