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

#ifndef MIMICCOMMLIB_TRANSPORT_H
#define MIMICCOMMLIB_TRANSPORT_H

#include "../message/Message.h"
#include "../serializer/Serializer.h"
#include "../DataTypes.h"



/**
 * Abstract class of transport layer
 *
 * TODO: Consider dividing into two abstract classes: MessageBasedTransport and Plain
 */
class Transport {

protected:
    /*
     * Serializer used for serialization/deserialization of messages
     */
    Serializer *serializer;

public:
    Transport(Serializer *serializer) : serializer(serializer) { }

    void setSerializer(Serializer *serializer) {
        Transport::serializer = serializer;
    }

    /**
     * Initialize server for receiving/sending data
     *
     * \param clients_number number of clients to get data from
     * \param replies map of client_id -> reply from the server
     */
    virtual void initServ(std::vector<std::string> paths) = 0;

    /**
     * Initialize client (prepare and connect to the server) for data interaction
     *
     * \param path path to working folder of the client
     * @todo consider merging into one init method
     */
    virtual void initClient(std::string path) = 0;

    /**
     * Establish connection to a specified address
     *
     * \param id id of the client to be connected
     */
    virtual int connectAddress (int id) = 0;

    /**
     * Accept a connection on a specified address
     *
     * \param id id of the client to accept connection from
     */
    virtual int acceptConnection(int id) = 0;

    /**
     * Disconnect client
     *
     * \param id id of the client to disconnect
     */
    virtual void closeConnection (int id) = 0;

    /**
     * Returns the address of the local server (remove?)
     */
    virtual char* getServerAddress () = 0;

    /**
     * Analogue to MPI_Probe return the size of the message waiting in the queue
     *
     * \param id id of the client
     * \param type DataType to be received
     */
    virtual int probe(int id, DataType type) = 0;

    /**
     * Destroy the port associated with the specified path
     */
    virtual void destroy(std::string path) = 0;

    /**
     * Send the data to a specified client
     *
     * \param data pointer to the array containing data to be sent
     * \param type type of the data to be sent
     * \param count number of elements in the data array
     * \param id of the receiving client
     */
    virtual void sendData(void *data, DataType type, int count, int id) = 0;

    /**
     * Receive data from the endpoint
     *
     * \param pointer to the array containing data to be sent
     * \param type type of the data to be received
     * \param count number of elements in the data array
     * \param id of the sending endpoint
     */
    virtual void receiveData(void *data_holder, DataType type, int count, int id) = 0;

};

#endif //MIMICCOMMLIB_TRANSPORT_H