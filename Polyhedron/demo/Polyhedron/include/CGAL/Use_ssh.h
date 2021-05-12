#ifndef USE_SSH_H
#define USE_SSH_H
#include <libssh/libsshpp.hpp>
#include <vector>

class QStringList;

namespace CGAL{
namespace ssh_internal{
//should be used inside a try/catch(ssh::SshException e)

//give an unitialized session.
bool establish_ssh_session(ssh_session& session,
                           const char *user,
                           const char *server,
                           const char *pub_key_path,
                           const char *priv_key_path,
                           const char *priv_key_password);

bool establish_ssh_session_from_agent(ssh_session& session,
                                      const char *user,
                                      const char *server,
                                      const char *pub_key_path);

void close_connection(ssh_session& session);

bool push_file(ssh_session& session,
               const char* dest_path,
               const char* filepath);
bool pull_file(ssh_session &session,
               const char *from_path,
               const char *to_path);

bool explore_the_galaxy(ssh_session &session,
                        QStringList &files);

}}
#endif
