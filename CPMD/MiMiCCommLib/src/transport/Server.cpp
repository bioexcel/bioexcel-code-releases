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

#include "Server.h"
#include "MPITransport.h"
#include <cstring>

int Server::init(std::vector<std::string> paths) {
    int clientsNumber = paths.size();
    std::vector<std::string> transformed_paths;
    for (int i = 0; i < clientsNumber; ++i) {
        transformed_paths.push_back(transform_path(paths[i]));
    }
    this->protocol->initServ(transformed_paths);
    this->clientsNumber = clientsNumber;

    for (int i = 0; i < clientsNumber; ++i) {
        client_list.push_back(new Client(i + 1, transformed_paths[i], protocol));
    }

    return 0;
}

void Server::disconnect(int dest) {
    protocol->closeConnection(dest);
}

void Server::destroy() {
    for (int i = 0; i < clientsNumber; ++i) {
        protocol->closeConnection(i + 1);
        protocol->destroy(client_list[i]->getPath());
    }
}

void Server::handshake() {
    for (int i = 0; i < clientsNumber; ++i) {
        protocol->acceptConnection(i + 1);
        int size = probe(i + 1, TYPE_CHAR);
        char *path = new char[size];
        protocol->receiveData(path, TYPE_CHAR, size, i + 1);
        std::string temp = path;
        temp = temp.substr(0, size);

        for (int j = 0; j < clientsNumber; j++) {
            if (std::strcmp(client_list[i]->getPath().c_str(), temp.c_str()) == 0) {
                int id = client_list[i]->getId();
                protocol->sendData(&id, TYPE_INT, 1, i + 1);
                break;
            }
        }
    }
}

int Server::send(void *data, int count, int destination, DataType type) {
    for (int i = 0; i < clientsNumber; i++) {
        if (client_list[i]->getId() == destination) {
            protocol->sendData(data, type, count, destination);
            break;
        }
    }
    return 0;
}

void Server::request(void *data, int count, int source, DataType type) {
    for (int i = 0; i < clientsNumber; i++) {
        if (client_list[i]->getId() == source) {
            protocol->receiveData(data, type, count, source);
            break;
        }
    }
}
