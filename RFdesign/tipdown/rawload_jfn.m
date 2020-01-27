function [dr] = rawload_jfn(pfile,frsize,nl,nframes,nslices,ncoils,frame,coils,slices,echoes)
% frsize -- number of data points acquire per TR
% nl     -- number of spiral leafs
% nframes -- number of time-frames
% nslices -- number of coils
% frame   -- get data for this time-frame only
% coils/slices/echoes -- get data for these coils/slices/echoes only
% $Id: rawload_jfn.m,v 1.1 2015/05/22 21:39:25 jfnielse Exp $


tothdrsize = 149788;   % MR750
rhbline = 1;

fip = fopen(pfile,'r');

% calculate size of data chunks
nviews = nl*nframes + rhbline;
echores = frsize*nviews;            % number of data points per echo (see pfilestruct.jpg)
sliceres = echores;        
coilres  = nslices * sliceres;      % number of data points per receive coil
totdatsize = ncoils*coilres*4;

% size of empty data
rhptsize = 2;  % data is stored in short int format (int16)
rhnecho = 1;
rhrawsize = 2*rhptsize*frsize*(nviews+0)*nslices*rhnecho;
rawoffset = 2*rhptsize*frsize*(0)*nslices*rhnecho;

% read data from file into array
% read only the desired coils,  slices, and echoes
dr = zeros(frsize, nl, length(slices), length(echoes), length(coils));
for coilind = 1 : length(coils)
	coil = coils(coilind);
	for sliceind = 1 : length(slices)
		slice = slices(sliceind);
		for echoind = 1 : length(echoes)
			echo = echoes(echoind);
			% offsetres = (coil-1)*coilres+ (slice-1)*sliceres + (echo-1)*echores + frsize*rhbline + (frame-1)*frsize*nl;
			offsetres = (coil-1)*coilres+ (slice-1)*sliceres + (echo-1)*echores + frsize*rhbline + (frame-1)*frsize*nl;
			offsetbytes = rawoffset + 2*rhptsize*offsetres;
			%if fseek(fip, -totdatsize+offsetbytes, 'eof') == -1
			if fseek(fip, tothdrsize+offsetbytes, 'bof') == -1
				error('Went beyond end of Pfile');
			end
			d = fread(fip, 2*frsize*nl, 'int16');
			%ftell(fip)
			d = complex(d(1:2:end), d(2:2:end));
			d = reshape(d, frsize, nl);
			dr(:, :, sliceind, echoind, coilind) = d;
		end;
	end;
end;

fseek(fip,0,'eof');
fprintf(1,'%s.m: Predicted header size = %d\n', mfilename, ftell(fip) - rhrawsize);
fclose(fip);

return;

% EOF
