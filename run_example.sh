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

# Download Slater-Koster files
SK_FILE_URL=https://dftb.org/fileadmin/DFTB/public/slako/auorg/auorg-1-1.tar.xz
if [ ! -d "auorg-1-1" ]; then
    curl --remote-name $SK_FILE_URL
    tar -xvf auorg-1-1.tar.xz
    rm auorg-1-1.tar.xz
fi

# Download database of SK parameters for TBMalt, if not downloaded already
DATABASE_PATH=example_dftb_parameters.h5
if [ ! -f $DATABASE_PATH ]; then
    python example_01_setup.py
fi

# Run both methods
python ase_dftb_qm9.py 1000
python tbmalt_run.py 1000 0 tbmalt_results.csv
python tbmalt_run.py 1000 1 tbmalt_results_individually.csv

# Make plot comparing them
python make_regression_table.py
Rscript qm9_regression.R
