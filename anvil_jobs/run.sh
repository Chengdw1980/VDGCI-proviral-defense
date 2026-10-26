cd /tmp/claude-0/-home-claude/1d89e018-43ce-5a67-b16c-37382686b1fc/scratchpad/disc48
while read a; do
  t=11; case $a in GCF_000027325.1|GCF_000027345.1) t=4;; esac
  [ -s faa/$a.faa ] || /tmp/env_df/bin/python /tmp/df5/callgenes.py g/$a.fna faa/$a.faa $t
  [ -d df/$a ] && [ -n "$(ls df/$a/*systems.tsv 2>/dev/null)" ] && continue
  PATH=/tmp/env_df/bin:$PATH /tmp/env_df/bin/defense-finder run --models-dir /tmp/dft/models --skip-model-version-check -w 2 -o df/$a faa/$a.faa > df/$a.log 2>&1 || echo "fail $a"
  echo -n "."
done < acc.txt
echo DONE
