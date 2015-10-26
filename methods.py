__author__="Vesselin Kolev <vlk@lcpe.uni-sofia.bg>"
__version__="20151001005601"


def get_ffnonbonded(cursor,ff_nonbonded):

   regexp='CREATE TABLE ff_nonbonded(aType TEXT, elNum INT, mass REAL, defCharge REAL, type TEXT, sigma REAL, eps REAL)'
   cursor.execute(regexp)

   regexp='CREATE INDEX ff_nonbonded_i ON ff_nonbonded(aType)'
   cursor.execute(regexp)

   flag=False

   f_obj=open(ff_nonbonded,"r")

   for line in f_obj:
      dummy=line.split()

      if flag:
         if len(dummy)>=1 and dummy[0]=='[':
            flag=False
      if len(dummy)>=3:
         if dummy[0]=='[' and dummy[1]=='atomtypes' and dummy[2]==']':
            flag=True
      if len(dummy)>=2 and flag:
         if (dummy[0][0] not in ["#",";","["]):
            regexp='INSERT INTO ff_nonbonded VALUES(?,?,?,?,?,?,?)'
            dummy[6]=dummy[6].split(';')[0]
            cursor.execute(regexp,(dummy[0],int(dummy[1]),float(dummy[2]),float(dummy[3]),dummy[4],float(dummy[5]),float(dummy[6])))

   f_obj.close()

   return True


def get_ffbonded(cursor,ff_bonded):

   regexp='CREATE TABLE ff_dihedrals(t1 TEXT, t2 TEXT, t3 TEXT, t4 TEXT, phi0 REAL, kd REAL, mult INT )'
   cursor.execute(regexp)

   regexp='CREATE INDEX ff_dihedrals_i ON ff_dihedrals(t1,t2,t3,t4)'
   cursor.execute(regexp)

   flag=False

   f_obj=open(ff_bonded,"r")

   for line in f_obj:
      dummy=line.split()

      if flag:
         if len(dummy)>=1 and dummy[0]=='[':
            flag=False
      if len(dummy)>=3:
         if dummy[0]=='[' and dummy[1]=='dihedraltypes' and dummy[2]==']':
            flag=True
      if len(dummy)>=2 and flag:
         if (dummy[0][0] not in ["#",";","["]):
            dummy[8]=dummy[8].split(';')[0]
            if int(dummy[4])==9:
               regexp='INSERT INTO ff_dihedrals VALUES(?,?,?,?,?,?,?)'
               cursor.execute(regexp,(dummy[0],dummy[1],dummy[2],dummy[3],float(dummy[5]),float(dummy[6]),int(dummy[7]),))

   f_obj.close()

   return True


def get_atoms_from_top(cursor,top_filename):

   regexp='CREATE TABLE atoms(aNum INT PRIMARY KEY, aType TEXT, resNum INT, resName TEXT, aName TEXT, cgrp INT,charge REAL, mass REAL, sigma REAL, eps REAL)'
   cursor.execute(regexp)

   f_obj=open(top_filename,"r")

   flag=False

   for line in f_obj:
      dummy=line.split()

      if flag:
         if len(dummy)>=1 and dummy[0]=='[':
            flag=False
      if len(dummy)>=3:
         if dummy[0]=='[' and dummy[1]=='atoms' and dummy[2]==']':
            flag=True
      if len(dummy)>=2 and flag:
         if (dummy[0][0] not in ["#",";","["]):
            regexp='SELECT * FROM ff_nonbonded WHERE aType=?'
            cursor.execute(regexp,(dummy[1],))
            sigma,eps=cursor.fetchall()[0][5:7]
            regexp='INSERT INTO atoms VALUES(?,?,?,?,?,?,?,?,?,?)'
            cursor.execute(regexp,(int(dummy[0]),dummy[1],int(dummy[2]),dummy[3],dummy[4],int(dummy[5]),float(dummy[6]),float(dummy[7]),sigma,eps,))

   f_obj.close()

   return True


def get_bonds_from_top(cursor,top_filename):

   regexp='CREATE TABLE bonds (left_a INT, right_a INT, left_a_type TEXT, right_a_type TEXT)'
   cursor.execute(regexp)

   f_obj=open(top_filename,"r")

   flag=False

   for line in f_obj:
      dummy=line.split()

      if flag:
         if len(dummy)>=1 and dummy[0]=='[':
            flag=False
      if len(dummy)>=3:
         if dummy[0]=='[' and dummy[1]=='bonds' and dummy[2]==']':
            flag=True
      if len(dummy)>=3 and flag:
         if (dummy[0][0] not in ["#",";","["]):
            if dummy[2].split(';')[0] is "1":
               regexp='INSERT INTO bonds(left_a,right_a) VALUES(?,?)'
               cursor.execute(regexp,(int(dummy[0]),int(dummy[1]),))

   f_obj.close()

   for i in ['left_a','right_a']:
      regexp='UPDATE bonds SET '+i+'_type=(SELECT aType FROM atoms WHERE aNum=bonds.'+i+')'
      cursor.execute(regexp)

   return True


def sqlite3_dihedral_resolver(cursor,aTypes):

   regexp='SELECT * FROM ff_dihedrals WHERE (t1=? and t2=? and t3=? and t4=?)'
   cursor.execute(regexp,(aTypes[0],aTypes[1],aTypes[2],aTypes[3],))

   result=cursor.fetchall()
   len_result=len(result)

   if len_result==0:
      regexp='SELECT * FROM ff_dihedrals WHERE (t1=? and t2=? and t3=? and t4=?)'
      cursor.execute(regexp,(aTypes[3],aTypes[1],aTypes[2],aTypes[0],))
      result=cursor.fetchall()
      len_result=len(result)

   if len_result==0:
      regexp='SELECT * FROM ff_dihedrals WHERE (t1=? and t2=? and t3=? and t4=?)'
      cursor.execute(regexp,("X",aTypes[1],aTypes[2],"X",))
      result=cursor.fetchall()
      len_result=len(result)

   if len_result==0:
      regexp='SELECT * FROM ff_dihedrals WHERE (t1=? and t2=? and t3=? and t4=?)'
      cursor.execute(regexp,("X",aTypes[2],aTypes[1],"X",))
      result=cursor.fetchall()
      len_result=len(result)

   if len_result==0:
      regexp='SELECT * FROM ff_dihedrals WHERE (t1=? and t2=? and t3=? and t4=?)'
      cursor.execute(regexp,(aTypes[0],aTypes[1],aTypes[2],"X",))
      result=cursor.fetchall()
      len_result=len(result)

   if len_result==0:
      regexp='SELECT * FROM ff_dihedrals WHERE (t1=? and t2=? and t3=? and t4=?)'
      cursor.execute(regexp,(aTypes[3],aTypes[1],aTypes[2],"X",))
      result=cursor.fetchall()
      len_result=len(result)

   if len_result==0:
      regexp='SELECT * FROM ff_dihedrals WHERE (t1=? and t2=? and t3=? and t4=?)'
      cursor.execute(regexp,(aTypes[3],aTypes[2],aTypes[1],"X",))
      result=cursor.fetchall()
      len_result=len(result)

   return len_result,result


def get_prop_dihedrals_from_top(cursor,top_filename):

   regexp='CREATE TABLE dihedrals(id INT, a1 INT, a2 INT, a3 INT, a4 INT, t1 TXT, t2 TXT, t3 TXT, t4 TXT, phi0 REAL, kd REAL, mult REAL)'
   cursor.execute(regexp)

   regexp='CREATE INDEX dihedrals_i ON dihedrals(t1,t2,t3,t4)'
   cursor.execute(regexp)

   f_obj=open(top_filename,"r")

   flag=False
   counter=1

   for line in f_obj:
      dummy=line.split()

      if flag:
         if len(dummy)>=1 and dummy[0]=='[':
            flag=False
      if len(dummy)>=3:
         if dummy[0]=='[' and dummy[1]=='dihedrals' and dummy[2]==']':
            flag=True
      if len(dummy)>=3 and flag:
         if (dummy[0][0] not in ["#",";","["]):
            if int(dummy[4])==9:
               regexp='INSERT INTO dihedrals(id,a1,a2,a3,a4) VALUES(?,?,?,?,?)'
               cursor.execute(regexp,(counter,int(dummy[0]),int(dummy[1]),int(dummy[2]),int(dummy[3]),))
               counter=counter+1

   f_obj.close()

   for i in xrange(4):
      regexp='UPDATE dihedrals SET t'+str(i+1)+'=(SELECT aType FROM atoms WHERE atoms.aNum=dihedrals.a'+str(i+1)+')'
      cursor.execute(regexp)

   return True


def get_coords_from_pdb(cursor,pdb_file):

   class Atom():
      pass

   f_obj=open(pdb_file,"r")

   for i in ['x','y','z']:
      regexp='ALTER TABLE atoms ADD COLUMN '+i+' REAL'
      cursor.execute(regexp)

   atom=Atom()

   for line in f_obj:
      line_spl=line.split()
      if line_spl[0]=='CRYST1':
         cryst1=line
      if line_spl[0]=='ATOM' or line_spl[0]=='HETATM':
         atom.ident=line[0:6]
         atom.num=int(line[6:11])
         atom.name=line[12:16]
         atom.altloc=line[16:17]
         atom.resname=line[17:20]
         atom.chainid=line[21:22]
         atom.ressec=line[22:26]
         atom.icode=line[26:27]
         atom.x=line[30:38]
         atom.y=line[38:46]
         atom.z=line[46:54]
         atom.occupancy=line[54:60]
         atom.tempfact=line[60:66]
         atom.segid=line[72:76]
         atom.elsym=line[76:78]
         atom.charge=line[78:80]

         regexp='UPDATE atoms SET x=?,y=?,z=? WHERE aNum=?'
         cursor.execute(regexp,(float(atom.x)/10.0,float(atom.y)/10.0,float(atom.z)/10.0,int(atom.num),))

   f_obj.close()

   return cryst1


def dump_atoms(cursor,atoms_file):

   f_obj=open(atoms_file,"w")

   regexp='SELECT COUNT(*) FROM atoms'
   cursor.execute(regexp)
   numbonds=int(cursor.fetchall()[0][0])

   line="%10d" % (int(numbonds))

   f_obj.write(line+'\n')

   regexp='SELECT * FROM atoms'
   cursor.execute(regexp)

   for i in cursor.fetchall():
#     202, u'O2', 20, u'ALA', u'OC1', 202, -0.8055, 16.0, 0.295992, 0.87864
      line="%10d %5s %10d %5s %5s %10d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f" %  (int(i[0]),str(i[1]),int(i[2]),str(i[3]),str(i[4]),int(i[5]),float(i[6]),float(i[7]),float(i[8]),float(i[9]),float(i[10]),float(i[11]),float(i[12]),)

      f_obj.write(line+'\n')

   f_obj.close()

   return True


def dump_bonds(cursor,bonds_file):

   f_obj=open(bonds_file,"w")

   regexp='SELECT COUNT(*) FROM bonds'
   cursor.execute(regexp)
   numbonds=int(cursor.fetchall()[0][0])

   line="%10d" % (int(numbonds))

   f_obj.write(line+'\n')

   regexp='SELECT * FROM bonds'
   cursor.execute(regexp)

   for i in cursor.fetchall():
      line="%10d %10d %5s %5s" %(int(i[0]),int(i[1]),str(i[2]),str(i[3]))
      f_obj.write(line+'\n')

   f_obj.close()

   return True


def dump_proper_dihedrals(cursor,prop_dih_file,prop_dih_par_file):

   f_obj=open(prop_dih_file,"w")

   regexp="SELECT COUNT(*) FROM dihedrals"
   cursor.execute(regexp)

   line="%10d" % int(cursor.fetchall()[0][0])

   f_obj.write(line+'\n')

   regexp='SELECT id,a1,a2,a3,a4,t1,t2,t3,t4 FROM dihedrals'
   cursor.execute(regexp)

   for i in cursor.fetchall():
      result=sqlite3_dihedral_resolver(cursor,i[5:10])
      line="%10d %10d %10d %10d %10d %2d " %(i[0],i[1],i[2],i[3],i[4],result[0])
      f_obj.write(line+'\n')
      if result[0]==0:
         print ':('

   f_obj.close()

   f_obj=open(prop_dih_par_file,"w")

   regexp='SELECT id,a1,a2,a3,a4,t1,t2,t3,t4 FROM dihedrals'
   cursor.execute(regexp)

   for i in cursor.fetchall():
      result=sqlite3_dihedral_resolver(cursor,i[5:10])
      if result[0]>0:
         for j in result[1]:
            line="%10d %12.6f %12.6f %2d" % (i[0],j[4],j[5]/2.0,j[6])
            f_obj.write(line+'\n')
      else:
         print ':('

   f_obj.close()

   return True
