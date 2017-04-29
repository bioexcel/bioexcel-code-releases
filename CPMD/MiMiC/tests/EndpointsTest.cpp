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


#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "../src/transport/Transport.h"
#include "../src/MCLMain.h"

using ::testing::Exactly;
using ::testing::_;
using ::testing::Return;

class MockTransport : public Transport {
public:
    MockTransport() : Transport(NULL) {}

    MOCK_METHOD1(initServ,
                 void(std::vector<std::string> paths));
    MOCK_METHOD1(initClient,
                 void(std::string path));
    MOCK_METHOD1(connectAddress,
                 int(int id));
    MOCK_METHOD1(acceptConnection,
                 int(int id));
    MOCK_METHOD1(closeConnection,
                 void(int id));
    MOCK_METHOD0(getServerAddress,
                 char*());
    MOCK_METHOD2(probe,
                 int(int id, DataType type));
    MOCK_METHOD1(destroy,
                 void(std::string path));
    MOCK_METHOD4(sendData,
                 void(void *data, DataType type, int count, int id));
    MOCK_METHOD4(receiveData,
                 void(void *data_holder, DataType type, int count, int id));
};

ACTION_P(SetArg0, value) {
    strcpy(static_cast<char*>(arg0), value);
}

TEST(ServerInit, ProperCalls) {
    std::vector<std::string> testAddresses;
    testAddresses.push_back("test1/");
    testAddresses.push_back("test2/");
    testAddresses.push_back("test3/");
    testAddresses.push_back("test4/");

    MockTransport* protocol = new MockTransport();
    Server server(protocol);

    EXPECT_CALL(*protocol, initServ(testAddresses)).Times(Exactly(1));

    server.init(testAddresses);

    delete(protocol);
}

TEST(ServerHandshake, ProperCalls) {
    std::vector<std::string> testAddresses;
    testAddresses.push_back("test1/");
    testAddresses.push_back("test2/");
    testAddresses.push_back("test3/");
    testAddresses.push_back("test4/");

    MockTransport* protocol = new MockTransport();
    Server server(protocol);
    server.addClient(new Client(0, testAddresses[0], protocol));
    server.addClient(new Client(1, testAddresses[1], protocol));
    server.addClient(new Client(2, testAddresses[2], protocol));
    server.addClient(new Client(3, testAddresses[3], protocol));
    server.set_clients_number(4);

    EXPECT_CALL(*protocol, acceptConnection(1)).Times(Exactly(1));
    EXPECT_CALL(*protocol, acceptConnection(2)).Times(Exactly(1));
    EXPECT_CALL(*protocol, acceptConnection(3)).Times(Exactly(1));
    EXPECT_CALL(*protocol, acceptConnection(4)).Times(Exactly(1));

    EXPECT_CALL(*protocol, probe(1, TYPE_CHAR))
            .Times(Exactly(1)).WillOnce(Return(6));
    EXPECT_CALL(*protocol, probe(2, TYPE_CHAR))
            .Times(Exactly(1)).WillOnce(Return(6));
    EXPECT_CALL(*protocol, probe(3, TYPE_CHAR))
            .Times(Exactly(1)).WillOnce(Return(6));
    EXPECT_CALL(*protocol, probe(4, TYPE_CHAR))
            .Times(Exactly(1)).WillOnce(Return(6));

    EXPECT_CALL(*protocol, receiveData(_, TYPE_CHAR, _, 1))
            .Times(Exactly(1)).WillOnce(SetArg0((char *)testAddresses[0].c_str()));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_CHAR, _, 2))
            .Times(Exactly(1)).WillOnce(SetArg0((char *)testAddresses[1].c_str()));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_CHAR, _, 3))
            .Times(Exactly(1)).WillOnce(SetArg0((char *)testAddresses[2].c_str()));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_CHAR, _, 4))
            .Times(Exactly(1)).WillOnce(SetArg0((char *)testAddresses[3].c_str()));

    EXPECT_CALL(*protocol, sendData(_, TYPE_INT, 1, 1)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(_, TYPE_INT, 1, 2)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(_, TYPE_INT, 1, 3)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(_, TYPE_INT, 1, 4)).Times(Exactly(1));

    server.handshake();

    delete(protocol);
}

TEST(ClientInit, ProperCalls) {
    std::vector<std::string> testAddresses;
    testAddresses.push_back("test1/");

    MockTransport* protocol = new MockTransport();
    Client client(protocol);

    EXPECT_CALL(*protocol, initClient(testAddresses[0])).Times(Exactly(1));

    client.init(testAddresses);

    delete(protocol);
}

TEST(ClientHandshake, ProperCalls) {
    std::vector<std::string> testAddresses;
    testAddresses.push_back("test1/");

    MockTransport* protocol = new MockTransport();
    Client client(1, testAddresses[0], protocol);
    client.addClient(new Server(protocol));

    EXPECT_CALL(*protocol, connectAddress(0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(_, TYPE_CHAR, 6, 0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_INT, 1, 0)).Times(Exactly(1));

    client.handshake();

    delete(protocol);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}