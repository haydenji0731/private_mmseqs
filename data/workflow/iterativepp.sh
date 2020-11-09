#!/bin/sh -e
# Iterative sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;

#targetdb_aln <- aln result
#targetdb <- profiles
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERYDB="$1"
TARGETDB="$2"
TMP_PATH="$4"

STEP=0

# processing
while [ $STEP -lt $NUM_IT ]; do
  if [ $STEP -eq 0 ]; then
    "$MMSEQS" search "$QUERYDB" "$TARGETDB" "$TMP_PATH/aln_$STEP" "$TMP_PATH" ${SEARCH_PAR} \
      || fail "Slicesearch died"
    if notExists "$TMP_PATH/profile_$STEP.dbtype"; then
      "$MMSEQS" expand2profile "$QUERYDB" "$TARGETDB" "$TMP_PATH/aln_$STEP" "$2_aln" "$TMP_PATH/profile_$STEP"  ${EXPAND_PAR} \
        || fail "Expand2Profile died"
    fi
    "$MMSEQS" profile2consensus "$TARGETDB" "$2_consensus" ${CONSENSUS_PAR} \
      || fail "Profile2Consensus died"
    TARGETDB="$2_consensus"
    # expand and then result2profile -> deal with the parameters as wel
  fi
  if [ $STEP -ge 1 ]; then
    if notExists "$TMP_PATH/pref_$STEP.dbtype"; then
      PARAM="PREFILTER_PAR_$STEP"
      eval TMP="\$$PARAM"
      if [ $STEP -eq 1 ]; then
        $RUNNER "$MMSEQS" prefilter "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_$STEP" ${TMP} \
          || fail "Prefilter died"
      else
        $RUNNER "$MMSEQS" prefilter "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_tmp_$STEP" ${TMP} \
          || fail "Prefilter died"
      fi
    fi
    #
    if [ $STEP -ge 1 ]; then
      if notExists "$TMP_PATH/pref_$STEP.dbtype"; then
        STEPPREV=$((STEP-1))
        "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_$STEP" "$TMP_PATH/aln_$STEPPREV" "$TMP_PATH/pref_$STEP" $SUBTRACT_PAR \
          || fail "Subtract died"
        "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_$STEP"
      fi
    fi
    # call alignment modulevi
    if notExists "$TMP_PATH/aln_tmp_$STEP.db"; then
      PARAM="ALIGNMENT_PAR_$STEP"
      eval TMP="\$$PARAM"
      if [ $STEP -eq 1 ]; then
        $RUNNER "$MMSEQS" align "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" ${TMP} \
          || fail "Alignment died"
      else
        $RUNNER "$MMSEQS" align "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_tmp_$STEP" ${TMP} \
          || fail "Alignment died"
      fi
    fi
    # merge alignment dbs
    if [ $STEP -ge 1 ]; then
      if notExists "$TMP_PATH/aln_$STEP.dbtype"; then
        STEPPREV=$((STEP-1))
        "$MMSEQS" mergedbs "$QUERYDB" "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEPPREV" "$TMP_PATH/aln_tmp_$STEP" \
          || fail "Mergedbs died"
        "$MMSEQS" rmdb "$TMP_PATH/aln_$STEPPREV"
        "$MMSEQS" rmdb "$TMP_PATH/aln_tmp_$STEP"
      fi
    fi
    # expand alignment dbs and create profiles
    if [ $STEP -ne $((NUM_IT - 1)) ]; then
      if notExists "$TMP_PATH/profile_$STEP.dbtype"; then
        PARAM="EXPANDPROFILE_PAR_$STEP"
        eval TMP="\$$PARAM"
        # make sure that $2_aln contains a backtrace
        $RUNNER "$MMSEQS" expand2profile "$QUERYDB" "$TARGETDB" "$TMP_PATH/aln_$STEP" "$2_aln" "$TMP_PATH/profile_$STEP" ${TMP} \
          || fail "Expand2Profile died"
      fi
      else
        if notExists "$3.dbtype"; then
        "$MMSEQS" expandaln "$QUERYDB" "$TARGETDB" "$TMP_PATH/aln_$STEP" "$2_aln" "$3" ${EXPANDALN_PAR} \
          || fail "Expandaln died"
      fi
    fi
  fi
  QUERYDB="$TMP_PATH/profile_$STEP"
  STEP=$((STEP+1))
done

if [ -n "$REMOVE_TMP" ]; then
  STEP=0
  while [ "$STEP" -lt "$NUM_IT" ]; do
    "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP" ${VERBOSITY}
    "$MMSEQS" rmdb "${TMP_PATH/aln_$STEP}" ${VERBOSITY}
    "$MMSEQS" rmdb "${TMP_PATH}/profile_$STEP" ${VERBOSITY}
    STEP=$((STEP+1))
  done
  rm -f "$TMP_PATH/iterativepp.sh"
fi