# Only execute if the folder does not exist
if [ ! -d "qm9_xyz" ]; then
    # Download QM9
    mkdir qm9_xyz
    cd qm9_xyz
    curl --location --remote-name --remote-header-name https://figshare.com/ndownloader/files/3195389
    INPATH=dsgdb9nsd.xyz.tar.bz2
    # Unzip a tar file compressed with bz2
    bzip2 -d $INPATH
    NEWINPATH=${INPATH%.bz2}
    tar -xvf $NEWINPATH
    # Remove the tar file
    rm $NEWINPATH
    cd ..
fi
