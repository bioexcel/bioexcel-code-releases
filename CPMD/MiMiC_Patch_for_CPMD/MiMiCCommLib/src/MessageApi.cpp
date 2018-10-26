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

#include "MessageApi.h"
#include "transport/MPITransport.h"
#include "MCLMain.h"

/**
 * External API of the library
 * NOTE! ALL API CALLS ARE BLOCKING!!!
 */

int MCL_init_server(char *paths_string, char delimeter) {
    MCLMain::getInstance().setProtocol(new MPITransport(MPI_COMM_SELF));
    return MCLMain::getInstance().initServer(paths_string, delimeter);
}

void MCL_init_client(char *path) {
    MCLMain::getInstance().setProtocol(new MPITransport(MPI_COMM_SELF));
    MCLMain::getInstance().initClient(path);
}

void MCL_send(void *data, int count, int data_type, int destination) {
    MCLMain::getInstance().send(data, count, data_type, destination);
}

void MCL_receive(void *buffer, int count, int data_type, int source) {
    MCLMain::getInstance().receive(buffer, count, data_type, source);
}

void MCL_destroy() {
    MCLMain::getInstance().destroy();
}

