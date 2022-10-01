# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 13:10:00 2022

@author: lawashburn
"""

from tkinter import *
import time
from random import randint
import threading

def five_seconds():
    time.sleep(5)
    my_label.config(text='5 seconds is up')

def rando():
    random_label.config(text=f'Random Number: {randint(1,100)}')

root = Tk()
root.title('Hi'
           )

root.geometry('500x400')

my_label = Label(root, text='Hello')
my_label.pack(pady=20)

mybutton1= Button(root, text='5 seconds',command = threading.Thread(target=five_seconds).start())
mybutton1.pack()

mybutton2= Button(root, text='Pick random number', command=rando)
mybutton2.pack()

random_label = Label(root,text='')
random_label.pack()


root.mainloop()