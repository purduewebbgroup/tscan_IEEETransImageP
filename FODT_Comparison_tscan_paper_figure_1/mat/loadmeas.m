function data=loadmeas(foo,N)
%function data=loadmeas(foo,N)

foo1=[foo '.dat'];
fid=datOpen(foo1,'r');
dims=[N 2];
[status, data]=read_float_array(fid, 'meas', length(dims), dims);
fclose(fid);

return;
