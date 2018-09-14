fid=datOpen('meas.dat','r');
dims=[312 2];
[status, data]=read_float_array(fid, 'meas', length(dims), dims);
fclose(fid);

fid=datOpen('meas2.dat','w+');
write_float_array(fid, 'meas', data, length(dims), dims);
fclose(fid);

fid=datOpen('meas2.dat','r');
[status, data2]=read_float_array(fid, 'meas', length(dims), dims);
fclose(fid);
