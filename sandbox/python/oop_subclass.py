#!/usr/bin/env python
class SchoolMember:
    '''any school member'''
    def __init__(self,name,age):
	'''initialize a school member'''
	self.name = name
	self.age = age
	print("Initialized SchoolMember: {}".format(self.name))
    def tell(self):
	'''tell details'''
	print('Name:"{}" Age:"{}"'.format(self.name,self.age)),
class Teacher(SchoolMember):
    '''teacher'''
    def __init__(self,name,age,salary):
	SchoolMember.__init__(self,name,age)
	self.salary = salary
	print("Initialized Teacher: {}".format(self.name))
    def tell(self):
	SchoolMember.tell(self)
	print('Salary:"{:d}"'.format(self.salary))
class Student(SchoolMember):
    '''student'''
    def __init__(self,name,age,marks):
	SchoolMember.__init__(self,name,age)
	self.marks = marks
	print("Initialized Student: {}".format(self.name))
    def tell(self):
	SchoolMember.tell(self)
	print('Marks:"{:d}"'.format(self.marks))
t = Teacher('DaGuo',30,1000)
s = Student('Yunfei',10,60)
print
for i in [t,s]:
    i.tell()
