#usage: seq2hla_forAlex hla/parent/path /target/path/

mkdir -p ${2}/HLA/ClassI
mkdir -p ${2}/HLA/ClassII

for f in $(find $1 -name *ClassI.HLAgenotype4digits ) ; do 
	myID=`dirname $f |sed 's/.*\///' `
	/hpc/users/akersn01/scripts/hla_grabber.pl $f | cut -f1 | sed "s/'//"  > ${2}/HLA/ClassI/${myID}.hla ; 
				# remove conf.  # remove apostrophes, they denote something...
done

for f in $(find $1 -name *ClassII.HLAgenotype4digits ) ; do 
	myID=`dirname $f |sed 's/.*\///' `
        /hpc/users/akersn01/scripts/hla_grabber.pl $f | cut -f1 | sed "s/'//"  > ${2}/HLA/ClassII/${myID}.hla ; 
                                # remove conf.  # remove apostrophes, they denote something...
done
