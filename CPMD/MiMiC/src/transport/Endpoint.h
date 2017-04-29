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

#ifndef MIMICCOMMLIB_ENDPOINT_H
#define MIMICCOMMLIB_ENDPOINT_H


#include <iostream>
#include "Transport.h"
#include "../DataTypes.h"

/**
 * Endpoint class representing clients involved in data communications
 */
class Endpoint {

protected:

    /**
     * Transport protocol used to send data
     */
    Transport* protocol;

    /**
     * List of known client (needed for server, but could be useful for p2p transports
     */
    std::vector<Endpoint*> client_list;

    /**
     * Endpoint id
     */
    int id;

    /**
     * Network address of the endpoint
     */
    std::string address;

    /**
     * Local path associated with the current endpoint
     */
    std::string path;

public:

    Endpoint(Transport *protocol) : protocol(protocol) { }


    Endpoint(int id, const std::string &path, Transport *protocol) :
            id(id), path(path), protocol(protocol) { }

    /**
     * Initialize endpoint
     *
     * \param paths paths to other endpoints which this one is
     * interacting with.
     */
    virtual int init(std::vector<std::string> paths) = 0;

    /**
     * Handshake procedure, assigns ids to clients
     */
    virtual void handshake() {}

    /**
     * Send data to a specific client (transport layer should handle data packaging)
     *
     * \param data pointer to the data array
     * \param count number of data entities
     * \param destination id of the client to receive data
     * \param type of the data to send
     */
    virtual int send(void *data, int count, int destination, DataType type) = 0;

    /**
     * Receive raw data (array of primitive type) from a specific client
     *
     * \param data pointer to the buffer
     * \param count number of data entities
     * \param source id of the client to receive data
     * \param type of the data to send
     */
    virtual void request(void *data, int count, int source, DataType type) = 0;

    /**
     * Probe message queue for length of the pending message
     *
     * \param client id of the client
     * \param type of the data to check
     */
    int probe(int client, DataType type) {
        return protocol->probe(client, type);
    }

    /**
     * Disconnect from a specified client
     */
    virtual void disconnect(int dest) = 0;

    /**
     * Destroy the endpoint
     */
    virtual void destroy () {}

    /**
     * Add new client to the clients list
     * @param client endpoint to add
     */
    void addClient(Endpoint* client) {
        client_list.push_back(client);
    }

    std::string transform_path(std::string original_path) {
        char sep = '/';
        std::string transformed_path = std::string(original_path);
        if (original_path[original_path.length() - 1] != sep) {
            transformed_path += sep;
        }
        return transformed_path;
    }

    int getId() {
        return id;
    }

    std::string getAddress() {
        return address;
    }


    std::string getPath() {
        return path;
    }

    void setId(int id) {
        Endpoint::id = id;
    }

    void setAddress(std::string address) {
        Endpoint::address = address;
    }

    void setPath(std::string path) {
        Endpoint::path = path;
    }

    const std::vector<Endpoint *> &getClient_list() const {
        return client_list;
    }
};


#endif //MIMICCOMMLIB_ENDPOINT_H
