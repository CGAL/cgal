#include <QtCore/QCoreApplication>
#include <QtCore/QCommandLineParser>
#include <QtCore/QCommandLineOption>
#include <QNetworkInterface>
#include <QApplication>
#include <QMessageBox>
#include "Server_ws.h"

#include <iostream>


EchoServer::EchoServer(quint16 port) :
  QObject(),
  m_pWebSocketServer(new QWebSocketServer(QStringLiteral("Echo Server"),
                                          QWebSocketServer::NonSecureMode, this))
{
  if (m_pWebSocketServer->listen(QHostAddress::Any, port)) {
    connect(m_pWebSocketServer, &QWebSocketServer::newConnection,
            this, &EchoServer::onNewConnection);
    connect(m_pWebSocketServer, &QWebSocketServer::closed, this, &EchoServer::closed);
  }
  QHostAddress local_host("0.0.0.0");

  //to avoid printing 127.0.0.1. Not realy sure it won't ever print the external ipv4 though.
  const QHostAddress &localhost = QHostAddress(QHostAddress::LocalHost);
  for (const QHostAddress &address: QNetworkInterface::allAddresses()) {
    if (address.protocol() == QAbstractSocket::IPv4Protocol && address != localhost)
    {
      local_host= address;
      break;
    }
  }
  QMessageBox mb(QMessageBox::NoIcon, "WS Server",
                 tr("WebSockets Server started.\nEnter the following address in\nyour Network Preferences to be able to join it :\n"
                    "ws://%1:%2").arg(local_host.toString()).arg(port), QMessageBox::Ok);
  mb.setTextInteractionFlags(Qt::TextSelectableByMouse);
  mb.exec();
}

EchoServer::~EchoServer()
{
  m_pWebSocketServer->close();
  qDeleteAll(m_clients.begin(), m_clients.end());
}

void EchoServer::onNewConnection()
{
  QWebSocket *pSocket = m_pWebSocketServer->nextPendingConnection();

  connect(pSocket, &QWebSocket::textMessageReceived, this, &EchoServer::processTextMessage);
  connect(pSocket, &QWebSocket::binaryMessageReceived, this, &EchoServer::processBinaryMessage);
  connect(pSocket, &QWebSocket::disconnected, this, &EchoServer::socketDisconnected);

  m_clients << pSocket;
}

void EchoServer::processTextMessage(QString message)
{
  QWebSocket *pClient = qobject_cast<QWebSocket *>(sender());
  for(auto *client : m_clients) {
    if(client != pClient)
      client->sendTextMessage(message);
  }
}

void EchoServer::processBinaryMessage(QByteArray message)
{
  QWebSocket *pClient = qobject_cast<QWebSocket *>(sender());
  if (pClient) {
    pClient->sendBinaryMessage(message);
  }
}

void EchoServer::socketDisconnected()
{
  QWebSocket *pClient = qobject_cast<QWebSocket *>(sender());
  if (pClient) {
    m_clients.removeAll(pClient);
    pClient->deleteLater();
  }
}

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  QCommandLineParser parser;
  parser.setApplicationDescription("WS Server");
  parser.addHelpOption();
  QCommandLineOption portOption(QStringList() << "p" << "port",
                                QCoreApplication::translate("main", "Port for echoserver [default: 1234]."),
                                QCoreApplication::translate("main", "port"), QLatin1Literal("1234"));
  parser.addOption(portOption);
  parser.process(a);
  int port = parser.value(portOption).toInt();
  EchoServer *server = new EchoServer(port);
  QObject::connect(server, &EchoServer::closed, &a, &QCoreApplication::quit);

  return a.exec();
}

