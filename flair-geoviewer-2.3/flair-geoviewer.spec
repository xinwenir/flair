%define _prefix /usr/local
%define bindir %{_prefix}/bin
%define prgdir %{_prefix}/flair
%define is_mandrake %(test -e /etc/mandrake-release && echo 1 || echo 0)
%define is_suse %(test -e /etc/SuSE-release && echo 1 || echo 0)
%define is_fedora %(test -e /etc/fedora-release && echo 1 || echo 0)
%define debug_package %{nil}

# to correctly generate an rpm for other distributions
%define _source_filedigest_algorithm 0
%define _binary_filedigest_algorithm 0
%define _binary_payload w9.gzdio

Name:    flair-geoviewer
Version: 2.3
Release: 0epy3
Packager: <paola.sala@orange.fr>
Prefix:  %{_prefix}
Source:  http://www.fluka.org/%{name}/%{name}-%{version}-%{release}.tgz
URL:     http://www.fluka.org/flair
BuildRoot: %{_tmppath}/%{name}-buildroot
License: Free for non-commercial non weapon related use

Summary: flair geometry viewer
Group: Applications/Engineering
Requires: flair = %{version}-%{release}
Requires: python3-imaging-tk
Requires: tk
#BuildRequires: python-devel
#BuildRequires: tk-devel
#BuildRequires: tcl-devel
AutoReqProv : no

%description
geoviewer is a 2D/3D geometry editor (viewer and debugger) for FLUKA geometries.

%prep
rm -rf $RPM_BUILD_ROOT
%setup -q

%build
make -j TARGET=%_target_cpu

%install
rm -Rf $RPM_BUILD_ROOT
make install ROOT=$RPM_BUILD_ROOT
#make install DESTDIR=$RPM_BUILD_ROOT%{prgdir}

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,0755)
#%doc AUTHORS BUGS README LICENSE ChangeLog
%{prgdir}/AUTHORS.geoviewer
%{prgdir}/BUGS.geoviewer
%{prgdir}/README.geoviewer
%{prgdir}/LICENSE.geoviewer
%{prgdir}/ChangeLog.geoviewer
%{prgdir}/*.so
%{prgdir}/usrbin2dvh
%{prgdir}/fonts/*.tga

%post
python3 << EOF
import os.path
import sys
import tkinter
libs = ["libtk8.5.so", "libtcl8.5.so", "libtk8.6.so", "libtcl8.6.so"]
found = [False]*len(libs)
pid = os.getpid()
try: f = open("/proc/%d/maps"%(pid),"r")
except: sys.exit(0)
for line in f:
   for i,lib in enumerate(libs):
      if not found[i] and line.find(lib)>=0:
         fn = line.split()[-1]
         path = os.path.dirname(fn)
         name = os.path.basename(fn)
         if name != lib:
            os.system("cd %s; ln -sf %s %s"%(path,name,lib))
         found[i] = True
f.close()
EOF

%changelog
* Sun May 5  2024 Paola Sala <paola.sala@mi.infn.it>
- Version: 2.3.0epy3
  - tagging 

* Fri Sep 1 2023 Paola Sala <paola.sala@mi.infn.it>
- Version: 2.3.0dpy3
  - Further bug fixing 

* Fri Sep 9 2022 Paola Sala <paola.sala@mi.infn.it>
- Version: 2.3.0cpy3
  - Bug fixing 

* Fri Jun 17 2022 Paola Sala <paola.sala@mi.infn.it>
- Porting to python3

* Fri Jun 18 2021 Paola Sala <paola.sala@mi.infn.it>
- Fixed license

* Sun Jan 1 2012 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Complete changelog can be found on the flair package
