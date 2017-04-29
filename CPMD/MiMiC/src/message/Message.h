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

#ifndef MIMICCOMMLIB_MESSAGE_H
#define MIMICCOMMLIB_MESSAGE_H

#include "BaseMessageData.h"
#include "FloatArrayData.h"
#include "SystemData.h"

/**
 * A message class representing messages sent between processes.
 */
struct Message {
public:
    /**
     * Message data
     */
    BaseMessageData* data;

    /**
     * Id of the process that sent a message
     */
    int sender_id;

    int message_code;

    Message() { }

    static const int OK_CODE = 1;

    static const int FAIL_CODE = -1;
};



#endif //MIMICCOMMLIB_MESSAGE_H