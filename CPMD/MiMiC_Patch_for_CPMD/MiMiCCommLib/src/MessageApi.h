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

#ifndef MIMICCOMMLIB_MESSAGESERVICE_H
#define MIMICCOMMLIB_MESSAGESERVICE_H


/**
 * \file MessageApi.h
 * \brief External API of the library
 *
 * Provides the C API to call library functions
 * NOTE! ALL API CALLS ARE BLOCKING!!!
 */

/**
 * \fn int MCL_init_server(char *paths, char delimeter)
 * Customizable init function. Uses arbitrary delimeter char
 *
 * \param paths local paths of all clients (needed for addresses sharing)
 * \param delimeter character which is used to extract individual paths
 */
extern "C"
int MCL_init_server(char *paths, char delimeter);

/**
 * \fn void MCL_init_client(char *path)
 * Initialize client endpoint
 *
 * \param path string containing the path in the file system to this client
 */
extern "C"
void MCL_init_client(char *path);

/**
 * \fn void MCL_send(void *data, int count, int data_type, int destination)
 * Send data to specified client
 *
 * \param data pointer to the buffer with data
 * \param count number of entities to send
 * \param data_type type of data to send
 * \param destination id of the client to receive data
 */
extern "C"
void MCL_send(void *data, int count, int data_type, int destination);

/**
 * \fn void MCL_receive(void *buffer, int count, int data_type, int source)
 * Receive data from a specified client
 *
 * \param buffer buffer to store data
 * \param count number of entities to receive
 * \param data_type type of data to send
 * \param source id of the client which is sending data
 */
extern "C"
void MCL_receive(void *buffer, int count, int data_type, int source);

/**
 * \fn void MCL_destroy()
 * Destroy the endpoint
 */
extern "C"
void MCL_destroy();


#endif //MIMICCOMMLIB_MESSAGESERVICE_H
