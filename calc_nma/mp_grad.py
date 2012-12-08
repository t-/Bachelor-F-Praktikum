#!/usr/bin/env python
import numpy as np
import os
import sys
import scipy.io as scpio

FUNIT= np.float64(512.3963) # unit conversion of frequency in THz 
FUNIT_long=np.float64(5.123963*10**14.0) # unit conversion of frequency in Hz
#FUNIT= np.float64(968.2887) # unit conversion of frequency in THz 
#FUNIT_long=np.float64(9.682887*10**14.0) # unit conversion of frequency in Hz
CUNIT=np.float64(29979245800.0) # c in cm/s [units]

class system():
    def __init__(self,filename):
        self.geom=self.read_xyz(filename)
        #self.dx=0.00529 # displacement for numerical gradients 0.01 bohr in units of Angstrom (molpro default)
        self.dx=0.00529
        self.hesse=np.zeros(shape=(self.geom.atom_count*3,self.geom.atom_count*3))
    class xyz():
        def __init__(self,atoms):
            self.atom_count=float(atoms[0][0])
            self.title=atoms[1][0]
            self.names=[]
            self.masses=[]
            self.coords=np.zeros(shape=(self.atom_count,3))
            self.init_coords(atoms[2:])
        def x(self,n):
            return self.coords[n][0]
        def y(self,n):
            return self.coords[n][1]
        def z(self,n):
            return self.coords[n][2]                
        def init_coords(self,atoms):
            for i,n in enumerate(atoms):
                self.names.append(n[0])
                self.coords[i][0]=float(n[1])
                self.coords[i][1]=float(n[2])
                self.coords[i][2]=float(n[3])
        def set_masses(self):
            for i in self.names:
                if i=='H':
                    self.masses.append(np.float64(1.00794))
                if i=='C':
                    self.masses.append(np.float64(12.0110))
                if i=='O':
                    self.masses.append(np.float64(15.99940))
                if i=='N':
                    self.masses.append(np.float64(14.00670))
                if i=='P':
                    self.masses.append(np.float64(30.97))
                if i=='S':
                    self.masses.append(np.float64(32.06))
            #print np.sum(np.array(self.masses))
        def write_xyz(self,filename='out.xyz'):
            f=open(filename,'w')
            print >> f, self.atom_count
            if self.title!='':
                print >> f, self.title
            else:
                print >> f, 'mp_grad'
            for i,j in enumerate(self.names):
                print >>f, j,' ',self.x(i),' ',self.y(i),' ',self.z(i)
            f.close()
        def write_xyz2(self,filename,atom_count,title,names,coords):
            f=open(filename,'w')
            print >> f, atom_count
            if title!='':
                print >> f, title
            else:
                print >> f, 'mp_grad'
            for i,j in enumerate(names):
                print >>f, j,' ',coords[i,0],' ',coords[i,1],' ',coords[i,2]
            f.close()
        def return_xyz(self):
            f=[]
            f.append(self.atom_count)
            if self.title!='':
                f.append(self.title)
            else:
                f.append('mp_grad')            
            for i,j in enumerate(self.names):
                f.append(str(j)+' '+str(self.x(i))+' '+str(self.y(i))+' '+str(self.z(i))  )
            return f

        def write_mpjob(self,filename,atom_count,title,names,coords):
            f=open( 'S2D2/'+filename + '.com','w')
            
            print >>f, "***, CAFFEINE cartesian coordinates (XYZ format)\nbasis=STO-3\ngeomtyp=xyz\nsymmetry,nosym\ngeometry={"
            print >> f, atom_count
            if title!='':
                print >> f, title
            else:
                print >> f, 'mp_grad'
            for i,j in enumerate(names):
                print >>f, j,' ',coords[i,0],' ',coords[i,1],' ',coords[i,2]
            print >>f,"}\nrhf\nE(1)=energy\ntable,E\ndigits,30\nsave,'"+filename+".dat'"
            
            f.close() 
                    
        def verbose(self):
            print self.coords
            print self.x(0),self.y(0),self.z(0)
    
    def read_xyz(self,filename):
        file2=[]
        for i in open(filename):
            line=i.strip('\n').split(' ')
            for i in xrange(line.count('')): line.remove('')
            file2.append(line)
        new_xyz=self.xyz(file2)
        new_xyz.set_masses()
        
        return new_xyz

    def Force_genV(self,column):
        dco=np.copy(self.geom.coords)
        dco.shape=(1,3*len(self.geom.coords))
        dco=dco[0]
        os.system('echo cd col'+str(column)+'>> run_mp.sh')
        for i in range(column,len(dco),1): # hessian is symmetrical, only elements below and on the diagonal are calculated!
            #self.write_submit_script('cols/col'+str(column)+'/', 'column'+str(i), str(column)+'_'+str(i) )
            #os.system('echo qsub column'+str(i)+'.sh >> run_mp.sh')
            #f=open('cols/col'+str(column)+'/column'+str(i)+'.com','w')
            f1=open('cols/col'+str(column)+'/column'+str(i)+'_1.xyz','w')
            f2=open('cols/col'+str(column)+'/column'+str(i)+'_2.xyz','w')
            f3=open('cols/col'+str(column)+'/column'+str(i)+'_3.xyz','w')
            f4=open('cols/col'+str(column)+'/column'+str(i)+'_4.xyz','w')
            #print >>f,"***,hessian column"+str(column)+" (XYZ format)\n"
            #print >>f,"basis=6-31g*\ngdirect\ngeomtyp=xyz\nsymmetry,nosym\n"
            di=np.zeros(len(dco))
            di[column]=self.dx

            ncalc=1
            for j in range(i,i+1):
                try:
                    dj=np.zeros(len(dco))
                    dj[j]=self.dx
                    print >>f1 , self.print_coor_str(dco+di+dj,ncalc)
                    ncalc+=1
                    print >>f2 , self.print_coor_str(dco-di+dj,ncalc)
                    ncalc+=1
                    print >>f3 , self.print_coor_str(dco+di-dj,ncalc)
                    ncalc+=1
                    print >>f4 , self.print_coor_str(dco-di-dj,ncalc)
                    ncalc+=1
                except Exception: 
                    print 'note: col too small', i,j
            #for ko in range(1,5):
             #os.system('cd cols/col'+str(column)+'/; babel '+'column'+str(i)+'_'+str(ko)+'.xyz'+' column'+str(i)+'_'+str(ko)+'.pdb;  cd ../../')
             #os.system('rm cols/col'+str(column)+'/column'+str(i)+'_'+str(ko)+'.xyz')
            #print >>f, "table,E\n digits,30\n save,'Energies_"+str(i)+".dat'\n"
        #os.system('echo cd ../ >> run_mp.sh')

    def V(self,geom,filename):#create molpro jobfile for energy calculation
        f=open(filename+'.com','w')
        print >>f, "***, CAFFEINE cartesian coordinates (XYZ format)\nmemory,300,m\nbasis=6-31g*\ngeomtyp=xyz\nsymmetry,nosym\ngeometry={"
        for i in geom.return_xyz():
            print >>f, i   
        #print >>f,"}\nrhf\nE(1)=energy\ntable,E\ndigits,30,30\nsave,'"+filename+".dat'"
        f.close()

    def read_V(self,column): #reads back the output of molpro energy runs
        results=[]
        for k in range(column,3*len(self.geom.coords),10): # hessian is symmetrical, only elements below and on the diagonal are calculated!
            try:
                for i in open('cols/col'+str(column)+'/energies_'+str(k)+'.dat','r'):
                    t=i.strip('\n').split(' ')
                    for j in range(t.count('')): t.remove('')
                    if len(t)>0 and t[0]!='E':
                        results.append(np.float64(t))
                        #print t
            except Exception:
                print 'cols/col'+str(column)+'/energies_'+str(k)+'.dat not FOunD!11'
        
        if ( len(results)- (3*len(self.geom.coords)-column)*4 ) == 0 :
            return np.array(results)
        else:
            print 'fatal error in column ',column,'\none of your jobs did not finish correctly!!'
            sys.exit(0)

    def fill_hesse(self,column,t):
        for i in range(0,len(t),4):
            self.hesse[int(column+i/4.0),column] = (t[i]-t[i+2]-t[i+1]+t[i+3])/(4.0*self.dx*self.dx)#/3.57106481 # 3.57106481 converts from 1/ang^2 to 1/bohr^2
            #print column,int(column+i/4.0), self.hesse[int(column+i/4.0),column]#/3.57106481
            if column+i/4.0!=column:
                self.hesse[column,column+i/4.0] = self.hesse[int(column+i/4.0),column]

    def nma_movie(self,v,omega,name):
        vn=v
        for i in range(0,1000,5):
            ii=float(i)/80.0;
            n=np.copy(self.geom.coords)
            n=n+0.35*vn*np.sin(2*ii)
            ii_int='%03d' % int(ii*10)
            self.geom.write_xyz2('out'+str(ii_int)+'.xyz',self.geom.atom_count,self.geom.title,self.geom.names,n)
        os.system('babel out*.xyz '+name)
        os.system('rm out*.xyz')

    def displace_struct(self,v,disp,tag='mode_1'):
        n=np.copy(self.geom.coords)
        n=n+disp*v
        self.geom.write_xyz2(tag+'displaced_'+str(disp)+'.xyz',self.geom.atom_count,self.geom.title,self.geom.names,n)

    def ppov(self,bb,geom,freq):
        print 'printing atom1 coordinates',str(1.0*geom.x(0)+2.3961050510)+','+str(1.0*geom.y(0)+5.0265011787)+','+str(1.0*geom.z(0)+47.6557502747)
        os.system("mkdir -p ray")
        for l in range(0,len(bb)):
            f=open('pov'+str(l)+'.pov','w')
            #print >>f, "camera {look_at <0.0,7.0,0.0>\n location <14.0 , 7.0 , 0.0>\n }"
            bbl=bb[:,l]
            print >>f, "background { color rgb <1, 1, 1> }"
            print >>f, "camera {direction<0.0,0.0,  -2.335>\n location <0.0 , 0.0 , 0.0>\n right 1.3333333731*x up y}"
            print >>f,"#default { finish{phong   -1.000 ambient    0.500 diffuse    0.450 phong_size 13.750000}}"
            print >>f, "light_source{<4000.0000,4000.0000,9999.0000>  rgb<1.0,1.0,1.0>}"
            
            for i in range(0,len(bb),3):
                a=0
                for j in range(0,3):
                    a+=bbl[i+j]*bbl[i+j]
                coor = str(1.0*geom.x(i/3)+0.08528781398)+','+str(1.0*geom.y(i/3)+0.35834771688)+','+str(1.0*geom.z(i/3)-48.3905013589)
                #print coor
                print >>f, "sphere{<"+coor+">, "+str(a)
                print >>f, "pigment{color rgb<1.0,0.0,0.0>}}"
            f.close()
            a=np.sqrt(np.abs(freq[l]))*FUNIT_long/CUNIT/(2.0*np.pi)
            os.system("cat stick.pov >>" +'pov'+str(l)+'.pov')
            #os.system("povray -D -W1920 -H1080 -I"+'pov'+str(l)+'.pov -Opov'+str(l)+'_'+str(a)+'.png')
            nnn='%.0f'%(a)
            #if l in [18,19]:
            os.system("povray -D -W1920 -H1200 -I"+'pov'+str(l)+'.pov -Opov'+str(l)+'_'+nnn+'.png')
            os.system("rm pov"+str(l)+".pov")
        
        
    def print_ev(self,freq,type2='wavenumber',printing='no'):
        if type2!='wavenumber':
            for number in freq:
                print np.sqrt(number)*FUNIT, 'THz'
                #return np.sqrt(number)*FUNIT
        else:
            for i,number in enumerate(freq):
                if number<0: print 'Im.',
                else: print 'Re.',
                if printing!='no':
                    print 'freq %4d % 6.2f THz, % 6.2f cm^-1' %(i+1,np.sqrt(np.abs(number))*FUNIT,np.sqrt(np.abs(number))*FUNIT_long/CUNIT/(2.0*np.pi))
                #return np.sqrt(np.abs(number))*FUNIT_long/CUNIT/(2.0*np.pi)

    def twonma(self,v1,v2,x1max,x2max,gmass,foldername):
        os.system('mkdir '+foldername)
        #transform equilibrium structure to mass weighted coordinates
        A = np.sqrt(gmass) * self.geom.coords.reshape(3*len(self.geom.coords))
        xrange1=np.linspace(-x1max,x1max,71)
        islast=False
        for i in xrange(len(xrange1)):
            x1=xrange1[i]
            #translate along v1 in mass weighted coordinates
            B=(A+x1*v1)
            if i==len(xrange1)-1: #needed for chain job
                islast=True
            self.twonma_gen_row(B,v2,x2max,gmass,i,foldername,islast)

    def twonma_gen_row(self,struct,v2,x2max,gmass,i,foldername,islast=False):
        xrange2=np.linspace(-x2max,x2max,71)
        f=open(foldername+'/struc'+str(i)+'.com','w')
        print >>f,"***,hessian column"+str(i)+" (XYZ format)\n"
        print >>f,"basis=6-31g*\ngdirect\ngeomtyp=xyz\nsymmetry,nosym\n"
        for j in xrange(len(xrange2)):
            x2=xrange2[j]
            B=(struct+x2*v2)
            B=B/np.sqrt(gmass)
            print >>f , self.print_coor_str(B.reshape(len(self.geom.coords),3),j)
        print >>f, "table,E\n digits,30\n save,'Energies_"+str(i)+".dat'\n"
        
        print >>f, "table,DX\n digits,30\n save,'DX_"+str(i)+".dat'\n"
        print >>f, "table,DY\n digits,30\n save,'DY_"+str(i)+".dat'\n"
        print >>f, "table,DZ\n digits,30\n save,'DZ_"+str(i)+".dat'\n"
        if islast:
            self.write_submit_script(foldername+'/','struc',i,nexts=False)
        else:
            self.write_submit_script(foldername+'/','struc',i,nexts=i+1)

    def print_coor_str(self,dco,ncalc):
        y = np.array(["%.10f" % w for w in dco.reshape(dco.size)])
        #y = y.reshape(dco.shape)
        y=y.reshape(self.geom.coords.shape)
        #s='geometry={\n'+str(self.geom.atom_count)+'\ngradgradgrad\n'
        s=str(int(self.geom.atom_count))+'\ngradgradgrad\n'
        for nr,i in enumerate(y):
            s+= str(self.geom.names[nr])+' '+' '.join(map(str, i))
            s+='\n'
        #s+='}\nrhf\nE('+str(ncalc+1)+')=energy\nDX('+str(ncalc+1)+')=DMX\nDY('+str(ncalc+1)+')=DMY\nDZ('+str(ncalc+1)+')=DMZ\n\n' # stupid molpro arrays start at 1, not zero! Thus the +1 shift
        return s

    def write_submit_script(self,path,name,nr,nexts=False):
        ff=open(path+name+str(nr)+".sh",'w')
        print >> ff,"#!/bin/bash"
        print >> ff,"#$ -S /bin/bash"
        print >> ff,"#$ -pe openmp_fast 1"
        print >> ff,"#$ -l mf=500M"
        #print >> ff,"#$ -l ds=1"
        print >> ff,"#$ -q *"
        print >> ff,"#$ -cwd"
        print >> ff,"#$ -l h_rt=23:13:37"
        print >> ff,"#$ -N Q"+str(nr)+'_'+path.replace('/','_')
        print >> ff,"#$ -M tgraen@gwdg.de"
        print >> ff,"#$ -m abe"
        print >> ff,". /etc/profile.d/modules.sh"
        print >> ff,"module load sge"
        print >> ff,"unset  PE_HOSTFILE"
        print >> ff,"/usr/local/molpro/bin/molpro -n 1 -m 60M -d $TMPDIR "+name+str(nr)+".com"
        print >> ff,"rm Q"+str(nr)+"* "+name+str(nr)+".out " + name+str(nr)+".xml "+name+str(nr)+".sh"
        if nexts and nexts not in [24,48]:
            print >> ff,"qsub "+name+str(nexts)+".sh >> id"+name+str(nexts)+".dat"
        ff.close()

    def twonma_read(self,foldername):
        a=np.ones(shape=(71,71))
        
        for i in range(0,71):
            energies=[]
            for line in open(foldername+'/energies_'+str(i)+'.dat','r'):
                ll=line.strip('').strip('\n').strip('E').split(' ')
                for il in xrange(ll.count('')):ll.remove('')
                if len(ll)==1:
                    energies.append(ll)
            if len(energies)!=71:
                print 'ERROR - read unfinished 2D potential: energies',i
            
            for v2 in range(0,71):
                #print np.size(a[:,1]), np.size(a[1,:]),'[v2,i], ',v2,i
                a[i,v2]=np.float64(energies[v2])
        print 'fin'
        scpio.savemat('2dpot_'+foldername+'.mat', mdict={'a': a})
    
    def write_hesse_input_files(self):
        self.V(self.geom,'out_com')
        os.system('echo ""> run_mp.sh')
        for i in range(0,3*len(self.geom.coords) ):
            os.system('mkdir -p cols/col'+str(i))
            self.Force_genV(column=i)

    def read_hesse_files(self):
        for i in range(0,3*len(self.geom.coords) ):
            t=self.read_V(column=i)
            self.fill_hesse(column=i,t=t)

    def write_hesse(self,filename):
        f=open(filename,'w')
        a=''
        for col in xrange(len(self.hesse[0])):
            for row in xrange(len(self.hesse)):
                if row<=col:
                    a=a+' % -2.7f' % self.hesse[row,col]
            print >>f,a
            a=''
            
    def calc_mass_vector(self):
        gmass=np.ones(3*len(self.geom.masses))
        for iii in xrange(3*len(self.geom.masses)):
            gmass[iii]=np.float64(self.geom.masses[int(iii/3.0)])
        #print 'gmass:',gmass
        return gmass
        
    def read_molpro_hesse(self,filename):
        m=[]
        for i in open(filename,'r'):
            a=i.strip('\n').split(' ')
            for k in range(a.count('')): a.remove('')
            for k in xrange(len(a)):
                a[k]=np.float64(a[k])
            if len(a)>0:
                m.append(a)
        return m
    
    def write_spec(self,aa,bb,gmass,filename):
        fev= open(filename,'w')
        tot=0.0 # total mass 3 times larger than single molecule
        for mode in xrange(len(bb)):
            si=0.0        
            for i in range(0,len(bb[:,mode]),3):
                si += gmass[i]*(bb[i,mode]**2.0+bb[i+1,mode]**2.0+bb[i+2,mode]**2.0)
            tot+=si
            print >>fev, mode,si, np.sqrt(np.abs(aa[mode]))*FUNIT 

    def inertia(self,filename):
        #=======================================================================
        #
        #    Calculates the moment of inertia tensor
        #
        #=======================================================================
        
        def shift_COM(g):
            Mass=0.0
            com=np.array([0.,0.,0.])
            for i in xrange(len(g.masses)):
                Mass+=g.masses[i]
                com[0]+=g.masses[i]*g.x(i)
                com[1]+=g.masses[i]*g.y(i)
                com[2]+=g.masses[i]*g.z(i)
            com/=Mass
            print 'Center of mass',com
            for i in xrange(len(g.masses)):
                g.coords[i]=g.coords[i]-com
            return g
        def rotate_axes(g,bb):
            shift_COM(g)
            for i in xrange(len(g.masses)):
                v=np.array([0.,0.,0.],dtype='float64')
                v[0] = np.dot(bb[:,0],g.coords[i])
                v[1] = np.dot(bb[:,1],g.coords[i])
                v[2] = np.dot(bb[:,2],g.coords[i])
                g.coords[i]=v
            return g
        
        def I11(m,geom):
            I=np.float64(0.0)
            for i in xrange(len(m)): I+=m[i]*(geom.y(i)**2.0 * geom.z(i)**2.0)
            return I
        def I22(m,geom):
            I=np.float64(0.0)
            for i in xrange(len(m)): I+=m[i]*(geom.x(i)**2.0 * geom.z(i)**2.0)
            return I
        def I33(m,geom):
            I=np.float64(0.0)
            for i in xrange(len(m)): I+=m[i]*(geom.x(i)**2.0 * geom.y(i)**2.0)
            return I
        def I12(m,geom):
            I=np.float64(0.0)
            for i in xrange(len(m)): I+=m[i]*(geom.x(i) * geom.y(i))
            return -I
        def I13(m,geom):
            I=np.float64(0.0)
            for i in xrange(len(m)): I+=m[i]*(geom.x(i) * geom.z(i))
            return -I
        def I23(m,geom):
            I=np.float64(0.0)
            for i in xrange(len(m)): I+=m[i]*(geom.y(i) * geom.z(i))
            return -I
            
        g=self.read_xyz('h2o.xyz')
        g=shift_COM(g)
        m=g.masses
        I=np.array([[I11(m,g),I12(m,g),I13(m,g)],[I12(m,g),I22(m,g),I23(m,g)],[I13(m,g),I23(m,g),I33(m,g)]])
        aa,bb=np.linalg.eig(I)
        g=rotate_axes(g,bb)
        I=np.array([[I11(m,g),I12(m,g),I13(m,g)],[I12(m,g),I22(m,g),I23(m,g)],[I13(m,g),I23(m,g),I33(m,g)]])
        
        a11=I[0,0]
        a12=I[0,1]
        a13=I[0,2]
        a21=I[1,0]
        a22=I[1,1]
        a23=I[1,2]
        a31=I[2,0]
        a32=I[2,1]
        a33=I[2,2]
        
        b11=np.float128(a22*a33-a23*a32)
        b12=np.float128(a13*a32-a12*a33)
        b13=np.float128(a12*a23-a13*a22)
        b21=np.float128(a23*a31-a21*a33)
        b22=np.float128(a11*a33-a13*a31)
        b23=np.float128(a13*a21-a11*a23)
        b31=np.float128(a21*a32-a22*a31)
        b32=np.float128(a12*a31-a11*a32)
        b33=np.float128(a11*a22-a12*a21)
        d=b11*(b22*b33-b23*b32)-b12*(b21*b33-b31*b23)+b13*(b21*b32-b22*b31)
        print '000000000'
        print b11/d,b12/d,b13/d
        print b21/d,b22/d,b23/d
        print b31/d,b32/d,b33/d
        print '000000000'
        
        aa,bb=np.linalg.eig(I)
        print 'inertia eigenval.:',aa
        
        print 'inertia ev1:',bb[:,0]
        print 'inertia ev2:',bb[:,1]
        print 'inertia ev3:',bb[:,2]
        
        print '___ev 0 1__orth.', np.dot(bb[:,0],bb[:,1])
        print '___ev 1 2__orth.', np.dot(bb[:,1],bb[:,2])
        print '___ev 0 2__orth.', np.dot(bb[:,0],bb[:,2])
        
        #print 'inverse check:\n',np.linalg.inv(bt1)*I*bt1-np.linalg.pinv(bt1)*I*bt1
        
        #def D1-3 translation, D4-6 rotation axis
        D1=np.zeros(3*len(m))
        D2=np.zeros(3*len(m))
        D3=np.zeros(3*len(m))

        # gen translational vectors
        for i in xrange(len(m)):
            D1[3*i+0]=np.sqrt(m[i])
            D2[3*i+1]=np.sqrt(m[i])
            D3[3*i+2]=np.sqrt(m[i])
        
        # gen rotational vectors -- broken code from gaussian white paper        
        #D4=np.zeros(3*len(m))
        #D5=np.zeros(3*len(m))
        #D6=np.zeros(3*len(m))
        #for i in xrange(len(m)):
        #    for j in xrange(3):
        #        D4[3*i+j]= ( np.dot(g.coords[i],bb[1,:]) * bb[j,2] - np.dot(g.coords[i],bb[2,:]) * bb[j,1] ) / np.sqrt(m[i])
        #        D5[3*i+j]= ( np.dot(g.coords[i],bb[2,:]) * bb[j,0] - np.dot(g.coords[i],bb[0,:]) * bb[j,2] ) / np.sqrt(m[i])
        #        D6[3*i+j]= ( np.dot(g.coords[i],bb[0,:]) * bb[j,1] - np.dot(g.coords[i],bb[1,:]) * bb[j,0] ) / np.sqrt(m[i])
        
        #rotate around Inertia Axis 1
        x=bb[0,0];y=bb[1,0];z=bb[2,0];
        x=1.0;y=0.0;z=0.0;
        dr1=np.array([[0,z,-y],[-z,0,x],[y,-x,0]])
        d1=np.zeros(3*len(m))
        for i in xrange(len(m)):
            for j in xrange(0,3):
                tm=np.dot(dr1,g.coords[i])
                d1[3*i+j]=tm[j]*np.sqrt(m[i])
        d1=d1/ np.linalg.norm(d1)
        #rotate around Inertia Axis 2        
        x=bb[0,1];y=bb[1,1];z=bb[2,1];
        x=0.0;y=1.0;z=0.0;
        dr2=np.array([[0,z,-y],[-z,0,x],[y,-x,0]])
        d2=np.zeros(3*len(m))
        for i in xrange(len(m)):
            for j in xrange(0,3):
                tm=np.dot(dr2,g.coords[i])
                d2[3*i+j]=tm[j]*np.sqrt(m[i])
        d2=d2/ np.linalg.norm(d2)
        #rotate around Inertia Axis 3        
        x=bb[0,2];y=bb[1,2];z=bb[2,2];
        x=0.0;y=0.0;z=1.0;
        dr3=np.array([[0,z,-y],[-z,0,x],[y,-x,0]])
        d3=np.zeros(3*len(m))
        for i in xrange(len(m)):
            for j in xrange(0,3):
                tm=np.dot(dr3,g.coords[i])
                d3[3*i+j]=tm[j]*np.sqrt(m[i])
        d3=d3/ np.linalg.norm(d3)
        # total rotation is sum of dd=d1+d2+d3 the x,y,z components of this vector are our new - rotation vectors 
        dd=d1+d2+d3
        ddx=np.zeros(3*len(m))
        ddy=np.zeros(3*len(m))
        ddz=np.zeros(3*len(m))
        for i in xrange(len(m)):
            ddx[3*i+0]=dd[3*i+0]
            ddy[3*i+1]=dd[3*i+1]
            ddz[3*i+2]=dd[3*i+2]
        
        print 'Check if rot axis are orthogonal d1,d2',np.dot(ddx,ddy)
        print 'Check if rot axis are orthogonal d1,d3',np.dot(ddx,ddz)
        print 'Check if rot axis are orthogonal d2,d3',np.dot(ddy,ddz)
        
        #normalize vectors 
        D1 = D1/ np.linalg.norm(D1)
        D2 = D2/ np.linalg.norm(D2)
        D3 = D3/ np.linalg.norm(D3)
        D4 = ddx/ np.linalg.norm(ddx)
        D5 = ddy/ np.linalg.norm(ddy)
        D6 = ddz/ np.linalg.norm(ddz)
        
        #transform matrix D
        D=np.zeros(shape=[3*len(m),3*len(m)])
        D[0]=D1
        D[1]=D2
        D[2]=D3
        D[3]=D4
        D[4]=D5
        D[5]=D6
        
        #create remaining set of vectors
        for i in range(6,3*len(m)):
            D[i][i]=1.0
        
        #gram-schmidt the remaining vectors
        def gs(D):
            for i in range(6,3*len(m)):
                Di=D[i]
                for j in range(0,i):
                    Di-= np.dot(D[j],D[i])*D[j]
                D[i]=Di
                #print '%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f ' %( np.abs(np.dot(D[i],D[0])),np.abs(np.dot(D[i],D[1])),np.abs(np.dot(D[i],D[2])),np.abs(np.dot(D[i],D[3])),np.abs(np.dot(D[i],D[4])),np.abs(np.dot(D[i],D[5])) )
        
        #for i in range(6,3*len(m)):
        #    print '%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f ' %( np.abs(np.dot(D[i],D[0])),np.abs(np.dot(D[i],D[1])),np.abs(np.dot(D[i],D[2])),np.abs(np.dot(D[i],D[3])),np.abs(np.dot(D[i],D[4])),np.abs(np.dot(D[i],D[5])),np.abs(np.dot(D[i],D[6])),np.abs(np.dot(D[i],D[7])) )
        print '\nsepsepsep\n'
        gs(D)
        print '    d3    d4    d5    '
        print 'd0% 4.8f % 4.8f % 4.8f' % (np.dot(D[0],D[3]),np.dot(D[0],D[4]),np.dot(D[0],D[5])  )
        print 'd1% 4.8f % 4.8f % 4.8f' % (np.dot(D[1],D[3]),np.dot(D[1],D[4]),np.dot(D[1],D[5])  )
        print 'd2% 4.8f % 4.8f % 4.8f' % (np.dot(D[2],D[3]),np.dot(D[2],D[4]),np.dot(D[2],D[5])  )
        print 'd3            % 4.8f % 4.8f'  % (                  np.dot(D[3],D[4]),np.dot(D[3],D[5])  )
        print 'd4                        % 4.8f'   % (                                    np.dot(D[4],D[5])  )
        
        return D.T
        
def main():
    print 'hello world'
    s=system('h2o.xyz')
    Dt=s.inertia('h2o.xyz')
    
    s.write_hesse_input_files()
    '''
    s.read_hesse_files()
    
    gmass=s.calc_mass_vector()
    corrmass= np.array(1/(np.sqrt(np.mat(gmass).T)*np.sqrt(np.mat(gmass))))
    s.write_hesse('hesse.dat')
    
    #In case a molpro hessian is available, read it here.
    #mhesse= s.read_molpo_hesse('t.dat')
    #aa,bb=np.linalg.eig(mhesse*corrmass)
    
    #Calc Eigenvectors and Eigenvalues
    fmwc=s.hesse*corrmass
    
    fint=np.dot(np.linalg.pinv(Dt),np.dot(fmwc,Dt))
    
    print 'translation', fint[0,0],fint[1,1],fint[2,2]
    print 'rotation   ', fint[3,3],fint[4,4],fint[5,5]
    for i in xrange(100):
            print fint[i,i]
    
    aa,bb=np.linalg.eig(fint[6:,6:])
    
    #s.print_ev(aa,type2='wavenumber',printing='yes')
    
    #write NMA spectrum    
    s.write_spec(aa,bb,gmass,'spec.dat')
    
    M=np.zeros(shape=[3*len(s.geom.masses),3*len(s.geom.masses)])
    for i in xrange(len(s.geom.masses)):
        for j in xrange(0,3):
            M[3*i+j,3*i+j]=1./np.sqrt(s.geom.masses[i])
    
    #L=np.zeros(shape=[3*len(s.geom.masses),3*len(s.geom.masses)])
    L=fint
    for i in range(0,len(np.array(bb))):
        for j in range(0,len(np.array(bb))):
            L[6+i,6+j]=np.array(bb)[i,j]
    
    F=np.zeros( 3*len(s.geom.masses) )
    for i in range(6,len(np.array(aa))):
        F[i]=aa[i]
        
    #calculate normal mode vectors in cartesian coordinates
    #aa,bb=np.linalg.eig(fint)
    lcart= np.dot(M,np.dot(Dt,L))
    
    M=np.zeros(shape=[3*len(s.geom.masses),3*len(s.geom.masses)])
    for i in xrange(len(s.geom.masses)):
        for j in xrange(0,3):
            M[3*i+j,3*i+j]=np.sqrt(s.geom.masses[i])

    #normalize normal mode vectors in cartesian coordinates
    #np.dot(M,lcart)
    for i in xrange(len(lcart)):
        N=np.sum(lcart[:,i]*lcart[:,i])
        if N!=0: lcart[:,i]=lcart[:,i]/np.sqrt(N)
    aa,bb=np.linalg.eig(fmwc)
    M=np.zeros(shape=[3*len(s.geom.masses),3*len(s.geom.masses)])
    for i in xrange(len(s.geom.masses)):
        for j in xrange(0,3):
            M[3*i+j,3*i+j]=1./np.sqrt(s.geom.masses[i])
    
    lcart=np.dot(M,bb)
    for i in xrange(len(lcart)):
        N=np.sum(lcart[:,i]*lcart[:,i])
        if N!=0: lcart[:,i]=lcart[:,i]/np.sqrt(N)
    #print Normal Mode Vectors Localization on Molecule
    aa,bb=np.linalg.eig(fmwc)
    #s.ppov(bb,s.geom,aa)
    print gmass
    targets=np.array([18,45,46,55,59,70,79,80,85,100])
    targets_xmax=[1.0472,1.3948,1.3937,1.4872,1.5367,1.5983,1.9107,1.9228,2.0207,3.6509]
    print np.sqrt(aa[targets])*FUNIT_long/CUNIT/(2.0*np.pi)
    for i in aa[targets]:
        print 'angular freq %4d omega=% 6.2f THz, f=% 6.2f cm^-1' %(i+1,np.sqrt(np.abs(i))*FUNIT,np.sqrt(np.abs(i))*FUNIT_long/CUNIT/(2.0*np.pi))
    #1/2 omega^2 x^2
    s=system('out_com.xyz')
    A = np.sqrt(gmass) * s.geom.coords.reshape(3*len(s.geom.coords))
    A = (A + 0.3*bb[:,18])/np.sqrt(gmass)
    s.geom.coords = A.reshape(len(s.geom.coords),3)
    s.geom.write_xyz2('Mana.xyz',s.geom.atom_count,s.geom.title,s.geom.names,s.geom.coords)

    s=system('out_com.xyz')
    A = np.sqrt(gmass) * s.geom.coords.reshape(3*len(s.geom.coords))
    A = (A - 0.3*bb[:,18])/np.sqrt(gmass)
    s.geom.coords = A.reshape(len(s.geom.coords),3)
    s.geom.write_xyz2('Mana2.xyz',s.geom.atom_count,s.geom.title,s.geom.names,s.geom.coords)
    
    
    #s.displace_struct( (bb[:,18]/np.sqrt(gmass)).reshape(s.geom.atom_count,3),0.0 )

    #print Normal Mode Vector Movies in kartesian coordinates
    #for i in xrange(len(aa)):    
    #    s.nma_movie((bb[:,i]/np.sqrt(gmass)).reshape(s.geom.atom_count,3),np.sqrt(np.abs(aa[i])),'mov_'+str(np.sqrt(np.abs(aa[i])))+'_'+str(i)+'.pdb')
    
    #2 NMA displacement potential
    targets=np.array([18,45,46,55,59,70,79,80,85,100])
    targets_xmax=[1.0472,1.3948,1.3937,1.4872,1.5367,1.5983,1.9107,1.9228,2.0207,3.6509]
    
    #for i,tar in enumerate(targets):
    #    for j in xrange(i+1,len(targets)):
    #        print 'writing,', targets[i],targets[j] 
    #        s.twonma(bb[:,targets[i]],bb[:,targets[j]],targets_xmax[i],targets_xmax[j],gmass,'E'+str(targets[i])+'_'+str(targets[j]))
    #for i,tar in enumerate(targets):
    #    for dis in xrange(-21,20,1):
    #        disp=dis*targets_xmax[i]/41.0
    #        s.displace_struct( (bb[:,targets[i]]/np.sqrt(gmass)).reshape(s.geom.atom_count,3),disp,str(targets[i])+'_')
    #
    #Prepare .m files for Matlab post-processing
    #li=      ['E18_100','E18_45','E18_46','E18_55','E18_59','E18_70','E18_79','E18_80','E18_85','E45_100','E45_46','E45_55','E45_59','E45_70','E45_79','E45_80','E45_85','E46_100','E46_55','E46_59','E46_70','E46_79','E46_80','E46_85','E55_100','E55_59','E55_70','E55_79','E55_80','E55_85','E59_100','E59_70','E59_79','E59_80','E59_85','E70_100','E70_79','E70_80','E70_85','E79_100','E79_80','E79_85','E80_100','E80_85','E85_100']
    #for i in li:
    #    s.twonma_read(i)
    '''
if __name__ == '__main__':
    main()
