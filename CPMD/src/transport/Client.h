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

#ifndef MIMICCOMMLIB_CLIENT_H
#define MIMICCOMMLIB_CLIENT_H


#include "../message/Message.h"
#include "Endpoint.h"

/**
 * Endpoint representing a client (slave) side
 */
class Client : public Endpoint {

public:
    Client(Transport* protocol) : Endpoint(protocol) { }

    Client(int id, const std::string &path, Transport *protocol) : Endpoint(id, path, protocol) { }

    int send(void *data, int count, int destination, DataType type);

    void request(void *data, int count, int source, DataType type);

    void handshake();

    int init(std::vector<std::string> paths);

    void disconnect(int dest);

    void destroy();

private:
    /**
     * Ugly, just ugly - but could be useful if at certain point we create a p2p
     * solution
     */
    static const int clientsNumber = 1;

};


#endif //MIMICCOMMLIB_CLIENT_H
