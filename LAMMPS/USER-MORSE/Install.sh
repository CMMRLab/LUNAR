run() { echo "$@"; $@; }

if [ "$1" = "1" ]; then
  if [ -e orig ]; then
    echo "Package is already installed";
  else
    rm -rf orig; mkdir orig;
    for i in $(find active -name '*.cpp' -o -name '*.h'); do
      if [ -e ../$(basename $i) ]; then run cp -p ../$(basename $i) orig; fi;
      run cp -p $i ..;
    done;
  fi;
elif [ "$1" = "0" ]; then
  if [ -e orig ]; then
    for i in $(find active -name '*.cpp' -o -name '*.h'); do
      run rm -f ../$(basename $i);
    done;
    for i in orig/*; do
      if [ -e $i ]; then run cp -p $i ..; fi;
    done;
    rm -rf orig;
  else
    echo "Package is not installed"
  fi;
fi;

