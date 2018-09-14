function data=loadsnr(foo,N)

foo1=[foo '.dat'];
fid=datOpen(foo1,'r');
dims=[N 1];
[status, data]=read_float_array(fid, 'snr', length(dims), dims);
fclose(fid);

return;
