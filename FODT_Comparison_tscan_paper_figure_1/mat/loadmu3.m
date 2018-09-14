function data=loadmu1(foo,var, dim1, dim2, dim3)
%function data=loadmu1(foo,var, dim1, dim2, dim3)

foo1=[foo '.dat'];
fid=datOpen(foo1,'r');
%dims=[2 129 33 17];
dims=[6 dim1 dim2 dim3];
[status, data]=read_float_array(fid, var, length(dims), dims);
fclose(fid);

return;
