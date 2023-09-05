pro quat_true_six, ts, tf, data_dir 

thisDevice = !D.NAME
SET_PLOT, 'Z'
LOADCT, 5
TVLCT, red, green, blue, /GET
SET_PLOT, thisDevice

interval = 5670.0 

sj = strarr(2)
sj[0] = data_dir 

set_plot, 'Z'

!P.MULTI = [0,1,3]
!P.CHARSIZE = 1.2
!X.TICKFORMAT = '(I6)'
!Y.TICKFORMAT = '(F8.4)'

sj[1] = '/diff_true_quat.dat'
del = strjoin(sj)
openr, 1, del
S1 = fstat(1) 
rec_line = S1.SIZE/63.0 ;
print, rec_line
cal = dblarr(4, rec_line);  dblarr(col_no,row_no) 
readf, 1, cal
close, 1 

sj[1] = '/ckg_3.dat'
sig = strjoin(sj)
openr, 1, sig
S1 = fstat(1)
rec_line2 = S1.SIZE/105.0 ;
print, rec_line2
err = dblarr(7, rec_line2);  dblarr(col_no,row_no)
readf, 1, err
close, 1

plot, cal(0,*), cal(1,*), /xstyle, /ystyle, $
		  title='TRUE - ESTI(LRS) on LRS x-axis', $
		  background = 255, color = 0, $
		  symsize = 0.1, psym = 2, $
                  xrange=[ts, tf], $
                  yrange=[-5, 5], $
                  ytitle='arcsec'

oplot, err(0,*), err(1,*), color = 120 ;
oplot, err(0,*), -err(1,*), color = 120 ;

plot, cal(0,*), cal(2,*),  /xstyle, /ystyle, $
		  title='TRUE - ESTI(LRS) on LRS y-axis', $
		  background = 255, color = 0, $
		  symsize = 0.1, psym = 2, $
                  xrange=[ts, tf], $
                  yrange=[-5, 5], $
                  ytitle='arcsec'

oplot, err(0,*), err(2,*), color = 120 ;
oplot, err(0,*), -err(2,*), color = 120 ;

plot, cal(0,*), cal(3,*), /xstyle, /ystyle, $
		  title='TRUE - ESTI(LRS) on LRS z-axis', $
		  background = 255, color = 0, $
		  symsize = 0.1, psym = 2, $
                  xrange=[ts, tf], $
                  yrange=[-10, 10], $
                  xtitle='second', $
                  ytitle='arcsec'

oplot, err(0,*), err(3,*), color = 120 ;
oplot, err(0,*), -err(3,*), color = 120 ;

; sj[1] = '/quat_true_six_' + ts + '_' + tf + '.png'
sj[1] = '/quat_true_six.png'
qsig = strjoin(sj)

image = tvrd()
thisImage = BYTSCL(image)
s = SIZE(thisImage)
image3d = BYTARR(3, s(1), s(2))
image3d(0, *, *) = red(thisImage)
image3d(1, *, *) = green(thisImage)
image3d(2, *, *) = blue(thisImage)
write_png, qsig, image3d

sj[1] = '/quat_true_six.sta'
stat = strjoin(sj)
openw, 1, stat 

result1 = moment(cal(1,*))
result2 = moment(cal(2,*))
result3 = moment(cal(3,*))

printf, 1, '             '
printf, 1, '                     X                 Y                 Z'

printf, 1, ' MAX error: ', max(abs(cal(1,*)),m1), $
                           max(abs(cal(2,*)),m2), $
                           max(abs(cal(3,*)),m3)
printf, 1, '              ', m1, '    ', m2, '    ', m3
printf, 1, ' MIN error: ', min(abs(cal(1,*)),n1), $
                           min(abs(cal(2,*)),n2), $
                           min(abs(cal(3,*)),n3)
printf, 1, '              ', n1, '    ', n2, '    ', n3

printf, 1, ' mean     : ', result1(0),  result2(0), result3(0)
printf, 1, ' std_dev  : ', sqrt(result1(1)), sqrt(result2(1)), sqrt(result3(1)) 
printf, 1, ' skewness : ', result1(2), result2(2), result3(2) 
printf, 1, ' curtosis : ', result1(3), result2(3), result3(3)

rms1 = sqrt(result1(0)*result1(0)+sqrt(result1(1))*sqrt(result1(1)))
rms2 = sqrt(result2(0)*result2(0)+sqrt(result2(1))*sqrt(result2(1)))
rms3 = sqrt(result3(0)*result3(0)+sqrt(result3(1))*sqrt(result3(1)))
printf, 1, ' RMS      : ', rms1, rms2, rms3

close, 1

end

