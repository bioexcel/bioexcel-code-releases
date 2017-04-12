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
#include "transport/Transport.h"
#include "transport/Client.h"
#include "transport/Server.h"
#include "transport/MPITransport.h"
#include <sstream>
#include <iterator>

/**
 * External API of the library
 * NOTE! ALL API CALLS ARE BLOCKING!!!
 */

Endpoint* endpoint;

int MCL_init_server(char *paths_string, char delimeter) {
    std::string merged_paths = std::string(paths_string);
    std::stringstream ss(merged_paths);
    std::string path;
    std::vector<std::string> client_paths;
    while (std::getline(ss, path, delimeter)) {
        if (!path.empty()) {
            client_paths.push_back(path);
        }
    }

    Transport* protocol = new MPITransport(MPI_COMM_SELF);
    Server* server = new Server(protocol);
    server->setId(0);
    endpoint = server;
    server->init(client_paths);
    server->handshake();
    return 0;
}

void MCL_init_client(char *path) {
    std::cout << "start initialization" << "\n";
    Transport* protocol = new MPITransport(MPI_COMM_SELF);
    std::string client_path = std::string(path);
    std::vector<std::string> paths;
    paths.push_back(client_path);
    Client* client = new Client(protocol);
    endpoint = client;
    std::cout << "initializing client" << "\n";
    client->init(paths);
    std::cout << "starting handshake" << "\n";
    client->handshake();
    std::cout << "Received id: " << client->getId() << "\n";
}

void MCL_send(void *data, int count, int data_type, int destination) {
    int temp_type = data_type;
    DataType type = static_cast<DataType>(temp_type);
    endpoint->send(data, count, destination, type);
}

void MCL_receive(void *buffer, int count, int data_type, int source) {
    int temp_type = data_type;
    DataType type = static_cast<DataType>(temp_type);
    endpoint->request(buffer, count, source, type);
}

void MCL_destroy() {
    endpoint->destroy();
}

