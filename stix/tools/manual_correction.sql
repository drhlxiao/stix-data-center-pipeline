update PCF set PCF_WIDTH=16 where PCF_NAME="NIX00123";
update PCF set PCF_WIDTH=8 where PCF_NAME="ZZPAD008";
update PCF set PCF_WIDTH=16 where PCF_NAME='ZZPAD016';
update PCF set PCF_WIDTH=24 where PCF_NAME='ZZPAD024';
update PCF set PCF_WIDTH=32 where PCF_NAME='ZZPAD032';
drop table if exists IDB;
create table IDB(
creation_datetime  datetime,
version  varchar(64) not null
);
insert into IDB (creation_datetime, version) values (current_timestamp, '2.26.28');

