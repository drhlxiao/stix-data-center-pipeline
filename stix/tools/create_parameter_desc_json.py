from core import idb

idb = idb.idb()
sql = 'select PCF_NAME, PCF_DESCR from PCF'
pcf = idb.execute(sql)
print("var PCF_DESCR={")
for desc in pcf:
    print('{}:"{}",'.format(desc[0], desc[1]))
sql = 'select CPC_PNAME, CPC_DESCR from CPC'
cpc = idb.execute(sql)
for desc in cpc:
    print('{}:"{}",'.format(desc[0], desc[1]))

print("};")

sql = 'select SCOS_NAME, SW_DESCR from SW_PARA'
scos = idb.execute(sql)
print("var SCOS_DESCR={")
for desc in scos:
    print('{}:"{}",'.format(desc[0], desc[1]))
print("};")
