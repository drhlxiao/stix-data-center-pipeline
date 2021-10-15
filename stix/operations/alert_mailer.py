#!/usr/bin/python3
# -*- coding: utf-8 -*-
""" This script is used to send stix instrument alerts (e.g. TM(5,2)) to STIX team

"""
import sys
import os
import smtplib
import socket
import datetime
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
import time
import requests

config = {
    'Email_from': 'stix_obs@fhnw.ch',
    'Email_user': '',
    'Email_pwd': '',
    'Email_server': 'lmailer.fhnw.ch',
    'port': 465
}

timeout = 60
socket.setdefaulttimeout(timeout)
DEBUG = False

HOST = 'https://pub023.cs.technik.fhnw.ch' if not DEBUG else 'http://localhost:5000'


def send_email(subject, content, conf, receivers):
    print('sending email...')
    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['From'] = conf['Email_from']
    msg['To'] = ", ".join(receivers)
    text = MIMEText(content, 'text')
    msg.attach(text)
    username = conf['Email_user']
    password = conf['Email_pwd']
    try:
        server = smtplib.SMTP_SSL(port=conf['port'])
        server.connect(conf['Email_server'])
    except socket.error as e:
        print('Failed to connect to the email server')
        return False
    except socket.timeout:
        print('server connection timeout')
        return False
    print('sending email...')
    if not DEBUG:
        server.sendmail(conf['Email_from'], receivers, msg.as_string())
    print("done!")
    server.quit()
    return True


def set_notification_sent(_id):
    print('updating notification...')
    url = '{}/request/operations/notifications/set-sent'.format(HOST)
    data = {'_id': _id}
    x = requests.post(url, data=data).json()
    print(x)


def main():
    url = '{}/request/operations/notifications'.format(HOST)
    notes = requests.get(url).json()
    for note in notes:
        _id = note['_id']
        content = note.get('content', None)
        title = note.get('title', None)
        receivers = note.get('user_emails', [])
        if content and receivers and title:
            if send_email(title, content, config, receivers):
                set_notification_sent(_id)
        else:
            print('ERROR: something is not invalid')
            print(note)
    print('done')


if __name__ == "__main__":
    main()
