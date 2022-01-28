from tkinter import Tk, Label, Button, Entry, IntVar, END, W, E, filedialog, StringVar, messagebox

import sys
import socket
import binascii

HOST = 'localhost'  # Standard loopback interface address (localhost)
PORTS = [9000, 9001]  # Port to listen on (non-privileged ports are > 1023)
SPEARATOR = '<-->'


def connect_socket(host, port):
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((host, port))
        return s
    except:
        return None


class MainFrame:
    def ask_filename(self, event):
        filename = filedialog.asksaveasfilename(
            initialdir=".",
            title="Set output filename",
            filetypes=(("raw data files", "*.dat"), ("all files", "*.*")))
        self.v.set(filename)
        filename = self.v.get()
        if self.f:
            self.f.close()
        self.f = open(filename, 'wb')

    def __init__(self, master):
        self.master = master
        master.title("MainFrame")
        self.f = None

        self.entered_number = 0

        self.msg = IntVar()
        self.msg.set('')
        self.total_label = Label(master, textvariable=self.msg)

        self.label = Label(master, text="Status:")
        self.label2 = Label(master, text="Filename:")

        self.v = StringVar()
        self.entry = Entry(master, textvariable=self.v)
        #self.entry.pack()
        self.v.set('Set output filename ...')
        self.entry.bind("<1>", self.ask_filename)

        self.start_button = Button(
            master, text="Start", command=lambda: self.start())
        self.stop_button = Button(
            master, text="Stop", command=lambda: self.stop())
        self.exit_button = Button(
            master, text="Exit", command=lambda: self.exit())

        # LAYOUT

        self.label.grid(row=0, column=0, sticky=W)
        self.total_label.grid(row=0, column=1, columnspan=2, sticky=E)

        self.entry.grid(row=1, column=1, columnspan=5, sticky=W)

        self.label2.grid(row=1, column=0, sticky=W)

        self.start_button.grid(row=2, column=0)
        self.stop_button.grid(row=2, column=1)
        self.exit_button.grid(row=2, column=2, sticky=W + E)
        self.running = True

    def stop(self):
        self.running = False
        if self.socket:
            self.socket.close()

    def exit(self):
        if self.f:
            self.close()
        MsgBox = messagebox.askquestion(
            'Exit Application',
            'Are you sure you want to exit the application',
            icon='warning')
        if MsgBox == 'yes':
            self.master.destroy()

    def start(self):
        for port in PORTS:
            self.socket = connect_socket(HOST, port)
            if self.socket:
                break
        if not self.socket:
            self.msg = 'Failed to connect to the socket...'
            return

        self.msg = 'socket connection established.'
        self.msg = 'waiting for packets ...'
        num = 0
        self.running = True
        self.buf = b''
        self.loop()

    def loop(self):
        if not self.running:
            return
        data = self.socket.recv(1024)
        if not data:
            return
        self.buf += data
        if self.buf.endswith(b'<-->'):
            data2 = self.buf.split()
            if self.buf[0:9] == 'TM_PACKET'.encode():
                data_hex = data2[-1][0:-4]
                data_binary = binascii.unhexlify(data_hex)
                self.f.write(data_binary)
                num += 1
                self.msg = 'Received TM ({}, {})(SPID {})  at {} '.format(
                    data2[3].decode(), data2[4].decode(), data2[1].decode(),
                    data2[2].decode())
            self.buf = b''
            self.master.after(0, daq)


root = Tk()
mainFrame = MainFrame(root)
root.mainloop()
