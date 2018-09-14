function data=loadphi(foo,var, K, dim1, dim2, dim3)
%function data=loadphi(foo,var, K, dim1, dim2, dim3)

foo1=[foo '.dat'];
fid=datOpen(foo1,'r');
dims=[K dim1 dim2 dim3 2];
[status, data]=read_float_array(fid, var, length(dims), dims);
fclose(fid);

return;

