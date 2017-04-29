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

#include "../message/Message.h"
#include "../message/FloatArrayData.h"
#include <ostream>
#include <istream>

#ifndef MIMICCOMMLIB_SERIALIZER_H
#define MIMICCOMMLIB_SERIALIZER_H

/**
 * Interface of serializer used to pack and unpack data from a parcel to be sent
 * over transport protocol
 */
class Serializer {

public:
    /**
     * Serialize data and put it into a stream
     *
     * \param msg message to be serialized
     * \param stream output strem in which to put data
     */
    virtual void serialize(Message msg, std::ostream *stream) = 0;

    /**
     * Deserialize data from stream
     *
     * \param stream with serialized data
     * \return Message contained in serialized data
     */
    virtual Message deserealize(std::istream *data) = 0;

};
#endif //MIMICCOMMLIB_SERIALIZER_H
