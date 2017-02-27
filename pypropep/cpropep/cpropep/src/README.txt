This is cpropep version 1.0, a program for determining the
theoretical performance of rocket propellant compositions.
Cpropep is released under the GPL and is written in ANSI C.
It has compiled successfully under the following platforms:

   . i686 Windows 98 (4.10.1998) running Microsoft Visual 
     C++ 6.0 Standard Edition with Microsoft Visual Studio
     Service Pack 2
   . i686 Red Hat Linux 6.2 (Kernel 2.2.14-12) with gcc
   . i586 Debian GNU/Linux 2.2 (Kernel 2.2.16) with gcc

To compile under linux, type make in the top-level directory.
Under MSVC++, create a new Win32 console project file, add 
all source and header files to the project, set the header
file directory to ..\lib\  (or copy all files from \lib to
\cpropep or change the #include statements in cpropep.c) and
build.


Run using 
    cpropep -h
for command-line help and usage.


Cpropep is Copyright (C) 2000 Antoine Lefebvre 
<antoinelefebvre@softhome.net>





