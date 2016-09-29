#! /usr/bin/awk -f

# ----------------------------------------------------------------

function module_filename(module_name) {
  return module_name".mod"
}

function find_in_array(name, arr, n_arr,   result)
{
  result = 0;

  if (n_arr > 0) {
    for (i = 1; i <= n_arr; i++) {
      if (name == arr[i]) {
        result = 1;
        break;
      }
    }
  }

  return result;
}

# ----------------------------------------------------------------

function file_init()
{
  # init for new file

  src_file = FILENAME
  obj_file = src_file
  sub(/\.[Ff]08$/, ".o", obj_file)
  file_defines_mod = 0

  FILES_f08[idx_src] = src_file; idx_src++
  FILES_obj[idx_obj] = obj_file; idx_obj++
}

function file_emit()
{
  # emit dependencies for previous file

  mod_file = module_filename(module)
  FILES_mod[idx_mod] = mod_file; idx_mod++
  printf "%s %s : %s", obj_file, mod_file, src_file
  for (i = 0; i < nd; i++) {
    print " \\"
    printf "  %s", module_filename(depends_on[i])
  }
  print ""
  printf("\t%s\n", COMMAND)
  print ""
}

function file_depend(  mod_file,   i)
{
  # Accumulate dependencies for previous file

  if (file_defines_mod > 0) {
    mod_file = module_filename(module)
    FILES_mod[idx_mod] = mod_file; idx_mod++
  }
  asort(FILES_mod);
  depend_rec = obj_file " " mod_file " : " src_file
  for (i = 0; i < nd; i++) {
    if (find_in_array(depends_on[i], sys_mods, n_sys_mods) == 0) {
      depend_rec = depend_rec " \\\n"
      depend_rec = depend_rec "  " module_filename(depends_on[i])
    }
  }
  depend_rec = depend_rec "\n" "\t" COMMAND "\n\n"

  return depend_rec
}

# ----------------------------------------------------------------

function lists_emit()
{
  # print ""
  # print "F08_SRCS", ":="
  # for (i = 1; i < idx_src; i++)
  #   print "F08_SRCS", "+=", FILES_f08[i]

  print "#-- DO NOT EDIT, AUTO GENERATED --",
      strftime("%a %b %d %H:%M:%S %Z %Y", systime())
  print ""
  print "F08_MODS", ":="
  for (i = 1; i < idx_mod; i++)
    print "F08_MODS", "+=", FILES_mod[i]
}

# ----------------------------------------------------------------

BEGIN {
  IGNORECASE = 1
  # FS = "[ \t]+"
  FS = "[ \t,]+"

  EMIT_DEP_LIST = ! EMIT_MOD_LIST

  first_init = 1
  idx_src = 1
  idx_mod = 1
  idx_obj = 1

  COMMAND = "$(F08_C) $(F08_FLAGS) -c $<"
  dependencies = ""

  n_sys_mods = split(EXCLUDE_MODS, sys_mods);
}

END {
  # file_emit()
  dependencies = dependencies file_depend()

  if (EMIT_MOD_LIST) {
    lists_emit()
    print ""
  }
  else {
    # print ""
    print dependencies
  }
}

# ----------------------------------------------------------------

(FNR == 1) {
  if (first_init == 1) {
    first_init = 0

    file_init()

    nd = 0
  }
  else {
    # file_emit()
    dependencies = dependencies file_depend()

    file_init()

    nd = 0
  }
}

/^([ \t]*)module([ \t]+)([a-zA-Z0-9_]+)/ {
  file_defines_mod = 1
  module = $2
}

/^([ \t]*)use([ \t]+)([a-zA-Z0-9_]+)(,?)/ {
  depends_on[nd] = $3
  nd++
}
