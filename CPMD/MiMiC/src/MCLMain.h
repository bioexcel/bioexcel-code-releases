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

#ifndef MIMICCOMMLIB_MCLMAIN_H
#define MIMICCOMMLIB_MCLMAIN_H

#include "transport/Endpoint.h"
#include "transport/Client.h"
#include "transport/Server.h"
#include <sstream>
#include <iterator>

/**
 * Singleton main class of the MCL - mostly needed in order to simplify testing
 */
class MCLMain {

private:
    Endpoint* endpoint;
    Transport* protocol;

    MCLMain() {};
    void operator=(MCLMain const&);

public:
    /**
     * Retrieve or create the singleton instance
     * @return singleton instance
     */
    static MCLMain& getInstance() {
        static MCLMain instance;
        return instance;
    }

    void setProtocol(Transport* protocol) {
        this->protocol = protocol;
    }

    /**
     * Customizable init function. Uses arbitrary delimiter char
     *
     * \param paths local paths of all clients (needed for addresses sharing)
     * \param delimiter character which is used to extract individual paths
     * \return error code - NOT IMPLEMENTED YET
     */
    int initServer(char *paths_string, char delimeter) {
        std::string merged_paths = std::string(paths_string);
        std::stringstream ss(merged_paths);
        std::string path;
        std::vector<std::string> client_paths;
        while (std::getline(ss, path, delimeter)) {
            if (!path.empty()) {
                client_paths.push_back(path);
            }
        }

        Server* server = new Server(this->protocol);
        server->setId(0);
        this->endpoint = server;

        server->init(client_paths);
        server->handshake();
        return 0;
    }

    /**
     * Initialize client endpoint
     *
     * \param path string containing the path in the file system to this client
     */

    void initClient(char *path) {
        std::string client_path = std::string(path);
        std::vector<std::string> paths;
        paths.push_back(client_path);

        Client* client = new Client(this->protocol);
        this->endpoint = client;

        client->init(paths);
        client->handshake();
    }

    /**
     * Send data to specified client
     *
     * \param data pointer to the buffer with data
     * \param count number of entities to send
     * \param data_type type of data to send
     * \param destination id of the client to receive data
     */
    void send(void *data, int count, int data_type, int destination) {
        int temp_type = data_type;
        DataType type = static_cast<DataType>(temp_type);
        this->endpoint->send(data, count, destination, type);
    }

    /**
     * Receive data from a specified client
     *
     * \param buffer buffer to store data
     * \param count number of entities to receive
     * \param data_type type of data to send
     * \param source id of the client which is sending data
     */
    void receive(void *buffer, int count, int data_type, int source) {
        int temp_type = data_type;
        DataType type = static_cast<DataType>(temp_type);
        endpoint->request(buffer, count, source, type);
    }

    /**
     * Destroy the endpoint
     */
    void destroy() {
        endpoint->destroy();
    }

    Endpoint *getEndpoint() const {
        return endpoint;
    }

};


#endif //MIMICCOMMLIB_MCLMAIN_H
