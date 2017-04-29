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

#include "Client.h"
#include "MPITransport.h"
#include "Server.h"

int Client::init(std::vector<std::string> paths) {
    protocol->initClient(transform_path(paths[0]));
    client_list.push_back(new Server(0, "", protocol));
    setPath(transform_path(paths[0]));

    return 0;
}

void Client::disconnect(int dest) {
    protocol->closeConnection(dest);
}

void Client::destroy() {
    disconnect(0);
}

void Client::handshake() {
    protocol->connectAddress(0);
    protocol->sendData((char *) getPath().c_str(), TYPE_CHAR, getPath().length(), 0);

    int id;
    protocol->receiveData(&id, TYPE_INT, 1, 0);
    this->setId(id);
}

int Client::send(void *data, int count, int destination, DataType type) {
    for (int i = 0; i < clientsNumber; i++) {
        if (client_list[i]->getId() == destination) {
            protocol->sendData(data, type, count, destination);
            break;
        }
    }
    return 0;
}

void Client::request(void *data, int count, int source, DataType type) {
    for (int i = 0; i < clientsNumber; i++) {
        if (client_list[i]->getId() == source) {
            protocol->receiveData(data, type, count, source);
            break;
        }
    }
}
