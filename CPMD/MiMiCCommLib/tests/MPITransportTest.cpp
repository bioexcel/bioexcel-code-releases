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
#include "MPITransport.h"
#include <mpi.h>

TEST(InitServer, ) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::string> testAddresses;
    testAddresses.push_back("./test1/");
    testAddresses.push_back("./test2/");

    if (rank == 0) {
        struct stat info;
        for (uint i = 0; i < testAddresses.size(); ++i) {
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
        for (uint i = 0; i < testAddresses.size(); ++i) {
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

    float sendFloat[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    float recvFloat[5];

    float sendDouble[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    float recvDouble[5];

    MPITransport transport(MPI_COMM_SELF);

    if (rank == 0) {
        struct stat info;
        for (uint i = 0; i < testAddresses.size(); ++i) {
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

        transport.sendData(sendFloat, TYPE_FLOAT, 5, 1);
        transport.sendData(sendFloat, TYPE_FLOAT, 5, 2);

        transport.receiveData(recvFloat, TYPE_FLOAT, 5, 1);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendFloat[i], recvFloat[i]);
        }
        transport.receiveData(recvFloat, TYPE_FLOAT, 5, 2);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendFloat[i], recvFloat[i]);
        }

        transport.sendData(sendDouble, TYPE_DOUBLE, 5, 1);
        transport.sendData(sendDouble, TYPE_DOUBLE, 5, 2);

        transport.receiveData(recvDouble, TYPE_DOUBLE, 5, 1);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendDouble[i], recvDouble[i]);
        }
        transport.receiveData(recvDouble, TYPE_DOUBLE, 5, 2);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendDouble[i], sendDouble[i]);
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

        transport.receiveData(recvFloat, TYPE_FLOAT, 5, 0);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendFloat[i], recvFloat[i]);
        }
        transport.sendData(recvFloat, TYPE_FLOAT, 5, 0);

        transport.receiveData(recvDouble, TYPE_DOUBLE, 5, 0);
        for (int i = 0; i < 5; ++i) {
            ASSERT_EQ(sendDouble[i], recvDouble[i]);
        }
        transport.sendData(recvDouble, TYPE_DOUBLE, 5, 0);

        transport.closeConnection(0);
        transport.destroy(testAddresses[rank - 1]);
    }
}

TEST(Probe, ) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::string> testAddresses;
    testAddresses.push_back("./test1/");
    testAddresses.push_back("./test2/");

    int send_size[2] = {5, 3};

    MPITransport transport(MPI_COMM_SELF);

    if (rank == 0) {
        struct stat info;
        for (uint i = 0; i < testAddresses.size(); ++i) {
            if( stat( testAddresses[i].c_str(), &info ) != 0 ) {
                mkdir(testAddresses[i].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }
        }
        MPITransport transport(MPI_COMM_SELF);
        transport.initServ(testAddresses);
        MPI_Barrier(MPI_COMM_WORLD);
        transport.acceptConnection(1);
        transport.acceptConnection(2);

        int size = transport.probe(1, TYPE_INT);
        ASSERT_EQ(send_size[0], size);
        int *recvBuffer = new int[size];
        transport.receiveData(recvBuffer, TYPE_INT, size, 1);

        delete(recvBuffer);
        size = transport.probe(2, TYPE_INT);
        ASSERT_EQ(send_size[1], size);
        recvBuffer = new int[size];
        transport.receiveData(recvBuffer, TYPE_INT, size, 2);

        transport.closeConnection(1);
        transport.closeConnection(2);
    } else {
        MPI_Barrier(MPI_COMM_WORLD);
        int *send_buffer = new int[send_size[rank - 1]];
        transport.initClient(testAddresses[rank - 1]);

        transport.connectAddress(0);

        transport.sendData(send_buffer, TYPE_INT, send_size[rank - 1], 0);
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