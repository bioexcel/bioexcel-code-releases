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

#include "gtest/gtest.h"
#include "../src/transport/MPITransport.h"
#include <mpi.h>

TEST(InitServer, ) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::string> testAddresses;
    testAddresses.push_back("./test1/");
    testAddresses.push_back("./test2/");

    if (rank == 0) {
        struct stat info;
        for (int i = 0; i < testAddresses.size(); ++i) {
            if( stat( testAddresses[i].c_str(), &info ) != 0 ) {
                mkdir(testAddresses[i].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }
        }

        MPITransport transport(MPI_COMM_SELF);
        transport.initServ(testAddresses);
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        MPI_Barrier(MPI_COMM_WORLD);
        std::string port_file = testAddresses[rank - 1].append(".portname");
        FILE* port_address = fopen(port_file.c_str(), "r");
        ASSERT_TRUE(port_address != NULL);
        char mpi_port_name[MPI_MAX_PORT_NAME];
        fscanf(port_address, "%s", mpi_port_name);
        fclose(port_address);
        std::string port(mpi_port_name);
        ASSERT_TRUE(!port.empty());
        remove(port_file.c_str());
    }
}

TEST(Connection, ) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::string> testAddresses;
    testAddresses.push_back("./test1/");
    testAddresses.push_back("./test2/");

    MPITransport transport(MPI_COMM_SELF);

    if (rank == 0) {
        struct stat info;
        for (int i = 0; i < testAddresses.size(); ++i) {
            if( stat( testAddresses[i].c_str(), &info ) != 0 ) {
                mkdir(testAddresses[i].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }
        }
        MPITransport transport(MPI_COMM_SELF);
        transport.initServ(testAddresses);
        MPI_Barrier(MPI_COMM_WORLD);
        int result;
        result = transport.acceptConnection(1);
        ASSERT_EQ(MPI_SUCCESS, result);
        result = transport.acceptConnection(2);
        ASSERT_EQ(MPI_SUCCESS, result);
        transport.closeConnection(1);
        transport.closeConnection(2);
    } else {
        MPI_Barrier(MPI_COMM_WORLD);
        transport.initClient(testAddresses[rank - 1]);

        int result = transport.connectAddress(0);
        ASSERT_EQ(MPI_SUCCESS, result);
        transport.closeConnection(0);
        transport.destroy(testAddresses[rank - 1]);
    }
}

TEST(PingPong, ) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::string> testAddresses;
    testAddresses.push_back("./test1/");
    testAddresses.push_back("./test2/");

    int sendBuffer[5] = {1, 2, 3, 4, 5};
    int recvBuffer[5];

    MPITransport transport(MPI_COMM_SELF);

    if (rank == 0) {
        struct stat info;
        for (int i = 0; i < testAddresses.size(); ++i) {
            if( stat( testAddresses[i].c_str(), &info ) != 0 ) {
                mkdir(testAddresses[i].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }
        }
        MPITransport transport(MPI_COMM_SELF);
        transport.initServ(testAddresses);
        MPI_Barrier(MPI_COMM_WORLD);
        transport.acceptConnection(1);
        transport.acceptConnection(2);

        transport.sendData(sendBuffer, TYPE_INT, 5, 1);
        transport.sendData(sendBuffer, TYPE_INT, 5, 2);

        transport.receiveData(recvBuffer, TYPE_INT, 5, 1);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendBuffer[i], recvBuffer[i]);
        }
        transport.receiveData(recvBuffer, TYPE_INT, 5, 2);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendBuffer[i], recvBuffer[i]);
        }

        transport.closeConnection(1);
        transport.closeConnection(2);
    } else {
        MPI_Barrier(MPI_COMM_WORLD);
        transport.initClient(testAddresses[rank - 1]);

        transport.connectAddress(0);

        transport.receiveData(recvBuffer, TYPE_INT, 5, 0);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendBuffer[i], recvBuffer[i]);
        }
        transport.sendData(recvBuffer, TYPE_INT, 5, 0);
        transport.closeConnection(0);
        transport.destroy(testAddresses[rank - 1]);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    if (numProcs != 3) {
        std::cout << "The test runs with 3 processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}