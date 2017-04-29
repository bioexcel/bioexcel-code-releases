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
using ::testing::Invoke;
using ::testing::WithArgs;
using ::testing::SetArgPointee;

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
    MCLMain main = MCLMain::getInstance();
    MockTransport* protocol = new MockTransport();
    main.setProtocol(protocol);
    std::vector<std::string> testAddresses;
    testAddresses.push_back("test1/");
    testAddresses.push_back("test2/");
    testAddresses.push_back("test3/");
    testAddresses.push_back("test4/");

    EXPECT_CALL(*protocol, initServ(testAddresses)).Times(Exactly(1));

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

    char *testString = (char *) "test1;test2;test3;test4";
    main.initServer(testString, ';');
    ASSERT_EQ(main.getEndpoint()->getId(), 0);
    ASSERT_EQ(main.getEndpoint()->getClient_list().size(), 4);
    const std::vector<Endpoint *> &clients = main.getEndpoint()->getClient_list();
    for (int i = 0; i < clients.size(); ++i) {
        ASSERT_EQ(testAddresses[i], clients[i]->getPath());
        ASSERT_EQ(i + 1, clients[i]->getId());
    }
    delete(protocol);
}

TEST(ClientInit, ProperCalls) {
    MCLMain main = MCLMain::getInstance();
    MockTransport* protocol = new MockTransport();
    main.setProtocol(protocol);
    std::string testAddress = "test/";
    char* testString = (char *) "test";

    EXPECT_CALL(*protocol, initClient(testAddress)).Times(Exactly(1));
    EXPECT_CALL(*protocol, connectAddress(0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(_, TYPE_CHAR, testAddress.length(), 0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_INT, 1, 0)).Times(Exactly(1));

    main.initClient(testString);

    ASSERT_EQ(testAddress, main.getEndpoint()->getPath());
    ASSERT_EQ(1, main.getEndpoint()->getClient_list().size());
    ASSERT_EQ(0, main.getEndpoint()->getClient_list()[0]->getId());

    delete(protocol);
}

TEST(ClientSend, ProperCalls) {
    MCLMain main = MCLMain::getInstance();
    MockTransport* protocol = new MockTransport();
    main.setProtocol(protocol);
    char* testString = (char *) "test";

    EXPECT_CALL(*protocol, initClient(_)).Times(Exactly(1));
    EXPECT_CALL(*protocol, connectAddress(0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(_, TYPE_CHAR, _, 0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_INT, 1, 0)).Times(Exactly(1));

    main.initClient(testString);

    int testData[5] = {1, 2, 3, 4, 5};

    EXPECT_CALL(*protocol, sendData(testData, TYPE_INT, 5, 0)).Times(Exactly(1));

    main.send(testData, 5, TYPE_INT, 0);

    delete(protocol);
}

TEST(ClientReceive, ProperCalls) {
    MCLMain main = MCLMain::getInstance();
    MockTransport* protocol = new MockTransport();
    main.setProtocol(protocol);
    char* testString = (char *) "test";

    EXPECT_CALL(*protocol, initClient(_)).Times(Exactly(1));
    EXPECT_CALL(*protocol, connectAddress(0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(_, TYPE_CHAR, _, 0)).Times(Exactly(1));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_INT, 1, 0)).Times(Exactly(1));

    main.initClient(testString);

    int testData[5];

    EXPECT_CALL(*protocol, receiveData(testData, TYPE_INT, 5, 0)).Times(Exactly(1));

    main.receive(testData, 5, TYPE_INT, 0);

    delete(protocol);
}

TEST(ServerSend, ProperCalls) {
    MCLMain main = MCLMain::getInstance();
    MockTransport* protocol = new MockTransport();
    main.setProtocol(protocol);
    char *testString = (char *) "test1;test2;test3;test4";

    EXPECT_CALL(*protocol, initServ(_)).Times(Exactly(1));
    EXPECT_CALL(*protocol, acceptConnection(_)).Times(Exactly(4));
    EXPECT_CALL(*protocol, probe(_, TYPE_CHAR)).Times(Exactly(4));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_CHAR, _, _)).Times(Exactly(4));

    main.initServer(testString, ';');

    int testData[5] = {1, 2, 3, 4, 5};

    EXPECT_CALL(*protocol, sendData(testData, TYPE_INT, 5, 1)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(testData, TYPE_INT, 5, 2)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(testData, TYPE_INT, 5, 3)).Times(Exactly(1));
    EXPECT_CALL(*protocol, sendData(testData, TYPE_INT, 5, 4)).Times(Exactly(1));

    main.send(testData, 5, TYPE_INT, 1);
    main.send(testData, 5, TYPE_INT, 2);
    main.send(testData, 5, TYPE_INT, 3);
    main.send(testData, 5, TYPE_INT, 4);

    delete(protocol);
}

TEST(ServerReceive, ProperCalls) {
    MCLMain main = MCLMain::getInstance();
    MockTransport* protocol = new MockTransport();
    main.setProtocol(protocol);
    char *testString = (char *) "test1;test2;test3;test4";

    EXPECT_CALL(*protocol, initServ(_)).Times(Exactly(1));
    EXPECT_CALL(*protocol, acceptConnection(_)).Times(Exactly(4));
    EXPECT_CALL(*protocol, probe(_, TYPE_CHAR)).Times(Exactly(4));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_CHAR, _, _)).Times(Exactly(4));

    main.initServer(testString, ';');

    int testData[5];

    EXPECT_CALL(*protocol, receiveData(testData, TYPE_INT, 5, 1)).Times(Exactly(1));
    EXPECT_CALL(*protocol, receiveData(testData, TYPE_INT, 5, 2)).Times(Exactly(1));
    EXPECT_CALL(*protocol, receiveData(testData, TYPE_INT, 5, 3)).Times(Exactly(1));
    EXPECT_CALL(*protocol, receiveData(testData, TYPE_INT, 5, 4)).Times(Exactly(1));

    main.receive(testData, 5, TYPE_INT, 1);
    main.receive(testData, 5, TYPE_INT, 2);
    main.receive(testData, 5, TYPE_INT, 3);
    main.receive(testData, 5, TYPE_INT, 4);

    delete(protocol);
}

TEST(DestroyServer, ProperCalls) {
    MCLMain main = MCLMain::getInstance();
    MockTransport* protocol = new MockTransport();
    main.setProtocol(protocol);
    char *testString = (char *) "test1;test2;test3;test4";

    std::vector<std::string> testAddresses;
    testAddresses.push_back("test1/");
    testAddresses.push_back("test2/");
    testAddresses.push_back("test3/");
    testAddresses.push_back("test4/");

    EXPECT_CALL(*protocol, initServ(_)).Times(Exactly(1));
    EXPECT_CALL(*protocol, acceptConnection(_)).Times(Exactly(4));
    EXPECT_CALL(*protocol, probe(_, TYPE_CHAR)).Times(Exactly(4));
    EXPECT_CALL(*protocol, receiveData(_, TYPE_CHAR, _, _)).Times(Exactly(4));

    main.initServer(testString, ';');

    EXPECT_CALL(*protocol, closeConnection(1)).Times(Exactly(1));
    EXPECT_CALL(*protocol, closeConnection(2)).Times(Exactly(1));
    EXPECT_CALL(*protocol, closeConnection(3)).Times(Exactly(1));
    EXPECT_CALL(*protocol, closeConnection(4)).Times(Exactly(1));

    EXPECT_CALL(*protocol, destroy(testAddresses[0])).Times(Exactly(1));
    EXPECT_CALL(*protocol, destroy(testAddresses[1])).Times(Exactly(1));
    EXPECT_CALL(*protocol, destroy(testAddresses[2])).Times(Exactly(1));
    EXPECT_CALL(*protocol, destroy(testAddresses[3])).Times(Exactly(1));

    main.destroy();
    delete(protocol);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}