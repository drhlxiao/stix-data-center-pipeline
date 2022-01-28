#!/usr/bin/env python3
# @title        : dump_socket_packets 
# @description  : A tool to receive packets via socket and dump the packets to a file 
# @author       : Hualin Xiao
# @date         : Feb. 11, 2019
#

import sys
import socket
import binascii
HOST = 'localhost'  # Standard loopback interface address (localhost)
PORTS = [9000,9001] # Port to listen on (non-privileged ports are > 1023)
SPEARATOR='<-->'
def connect_socket(host,port):
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((host, port))
        return s
    except:
        return None

def run(filename):
    for port in PORTS:
        s = connect_socket(HOST, port)
        if s:
            break
    if not s:
        print('Failed to connect to the socket...')
        return 

    print('socket connection established.')
    print('waiting for packets ...')
    with open(filename,'wb') as f:
        num = 0
        try:
            buf=b''
            while True:
                data = s.recv(1024)
                if not data:
                    break
                buf += data
                if buf.endswith(b'<-->'):
                    data2=buf.split()
                    if buf[0:9] == 'TM_PACKET'.encode():
                        data_hex=data2[-1][0:-4]
                        data_binary = binascii.unhexlify(data_hex)
                        f.write(data_binary)
                        num += 1
                        print('Received: TM ({}, {})(SPID {})  at {} '.format(data2[3].decode(), 
                            data2[4].decode(), data2[1].decode(),data2[2].decode()))
                    buf=b''
        except KeyboardInterrupt:
            f.close()
            print('{} packets written to {}'.format(num, filename))

        finally:
            s.close()


if __name__=='__main__':
    if len(sys.argv) < 2:
        print('dump_socket_packets  <filename>')
    else:
        run(sys.argv[1])
