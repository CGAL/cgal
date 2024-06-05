#ifndef SERVER_WS_H
#define SERVER_WS_H

#include  <QObject>
#include <QWebSocketServer>
#include <QWebSocket>

class EchoServer : public QObject
{
  Q_OBJECT
public:
  explicit EchoServer(quint16 port);
  ~EchoServer();


Q_SIGNALS:
  void closed();

private Q_SLOTS:
  void onNewConnection();
  void processTextMessage(QString message);
  void processBinaryMessage(QByteArray message);
  void socketDisconnected();

private:
  QWebSocketServer *m_pWebSocketServer;
  QList<QWebSocket *> m_clients;
};

#endif // SERVER_WS_H
