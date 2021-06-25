new_files=$(find /data/gfts/solsoc/from_soc/ -name "solo_ANC*zip" -newermt "-24 hours" )
echo "new spice kernel files:"
echo $new_files
if [ -n $new_files ]; then
	echo "unzip spice kernels"
	for f in $new_files
	do 
		/usr/bin/unzip -o $f -d /data/gfts/data/spice/latest/
	done
	ln -s /data/gfts/solsoc/from_soc/solo_ANC*  /data/pub/data/spice/archive/ 2>>/dev/null
	sleep 1
	echo "restarting parser pipeline"
	kill `ps -eo pid,args --cols=10000 | awk '/daemon.py/ && $1 != PROCINFO["pid"] { print $1 }'`
	cd /opt/stix/parser
	nohup  python3 stix/pipeline/daemon.py&
	sleep 1
	cd /opt/stix/web
	echo "restart gunicorn"
	pkill gunicorn
	gunicorn --bind 127.0.0.1:8001 wsgi:app --daemon --threads=4 --workers=1 --reload --log-file ./out.log --worker-class=gthread

fi

