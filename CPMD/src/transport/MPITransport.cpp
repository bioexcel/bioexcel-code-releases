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


#include <stdexcept>
#include <cstdio>
#include "MPITransport.h"

std::string MPITransport::FILENAME = ".portname";

void MPITransport::initServ(std::vector<std::string> paths) {
    unsigned int client_number = (unsigned int) paths.size();
    port = new mpi_port_name[client_number + 1];
    intercomm = new MPI_Comm[client_number + 1];
    for (int i = 0; i < client_number; ++i) {
        MPI_Open_port(MPI_INFO_NULL, port[i + 1]);
        std::string filepath = paths[i].append(FILENAME);
        remove(filepath.c_str());
        FILE *port_file = fopen(filepath.c_str(), "w");
        fprintf(port_file, port[i + 1]);
        fclose(port_file);
    }
}

void MPITransport::initClient(std::string path) {
    bool file_exists = false;
    FILE *port_address = NULL;
    std::string port_file = path.append(FILENAME);
    port = new mpi_port_name[1];
    intercomm = new MPI_Comm[1];
    while (!file_exists) {
        port_address = fopen(port_file.c_str(), "r");
        if (port_address != NULL) {
            file_exists = true;
        }
    }
    fscanf(port_address, "%s", port[0]);
    fclose(port_address);
}

int MPITransport::connectAddress(int id) {
    MPI_Comm_connect(port[id], MPI_INFO_NULL, 0, host_comm, &intercomm[id]);
    return 0;
}

int MPITransport::acceptConnection(int id) {
    MPI_Comm_accept(port[id], MPI_INFO_NULL, 0, host_comm, &intercomm[id]);
    return 0;
}

void MPITransport::closeConnection(int id) {
    MPI_Comm_disconnect(&intercomm[id]);
}

char *MPITransport::getServerAddress() {
    return port[0];
}

void MPITransport::sendData(void *data, DataType type, int count, int id) {
    MPI_Ssend(data, count, pick_mpi_type(type), 0, 0, intercomm[id]);
}

void MPITransport::receiveData(void *data_holder, DataType type, int count, int id) {
    MPI_Status status;
    MPI_Recv(data_holder, count, pick_mpi_type(type), 0, MPI_ANY_TAG, intercomm[id], &status);
}

int MPITransport::probe(int id, DataType type) {
    int size;
    MPI_Status status;
    MPI_Probe(0, MPI_ANY_TAG, intercomm[id], &status);
    MPI_Get_count(&status, pick_mpi_type(type), &size);
    return size;
}

void MPITransport::destroy(std::string path) {
    std::string filepath = path.append(FILENAME);
    remove(filepath.c_str());
}

MPI_Datatype MPITransport::pick_mpi_type(DataType type) {
    MPI_Datatype send_type;
    switch (type) {
        case TYPE_CHAR:
            send_type = MPI_CHARACTER;
            break;

        case TYPE_INT:
            send_type = MPI_INT;
            break;

        case TYPE_DOUBLE:
            send_type = MPI_DOUBLE;
            break;

        case TYPE_FLOAT:
            send_type = MPI_FLOAT;
            break;

        default:
            send_type = MPI_CHAR;
    }
    return send_type;
}
