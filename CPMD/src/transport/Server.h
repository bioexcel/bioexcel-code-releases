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

#ifndef MIMICCOMMLIB_SERVER_H
#define MIMICCOMMLIB_SERVER_H


#include "Endpoint.h"
#include "Client.h"

/**
 * Endpoint representing server(master)
 */
class Server : public Endpoint {

private:
    int clientsNumber;

public:
    Server(Transport *protocol) : Endpoint(protocol) { }

    Server(int id, const std::string &path, Transport *protocol) : Endpoint(id, path, protocol) { }

    int init(std::vector<std::string> paths);

    void disconnect(int dest);

    void destroy();

    void handshake();

    int send(void *data, int count, int destination, DataType type);

    void request(void *data, int count, int source, DataType type);

    void set_clients_number(int clients_number) {
        Server::clientsNumber = clients_number;
    }

    int getClientsNumber() {
        return clientsNumber;
    }
};


#endif //MIMICCOMMLIB_SERVER_H
