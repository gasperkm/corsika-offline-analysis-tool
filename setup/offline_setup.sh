# Environmental variables for ADST reader 

# Old/stable version of Offline/APE
olddir=/install/path/for/old/ape/version
# New/trunk version of Offline/APE
newdir=/install/path/for/new/ape/version

if [ "$1" == "old" ]; then
   cd $olddir
   ./aperc_setup.sh
   echo "To setup environmental variables, use:"
   echo "   eval \`ape sh offline\`"
elif [ "$1" == "new" ]; then
   cd $newdir
   ./offline_devel.sh aperc
   echo "To setup environmental variables, use:"
   echo "   source /data0/gkukec/Programi/offline-trunk/env_offline.sh"
else
   echo "No arguments supplied (old/new)."
fi
