pro quat_true_quest, ts, tf, data_dir 

interval = 5670.0 

sj = strarr(2)
sj[0] = data_dir 

set_plot, 'Z'

!P.MULTI = [0,1,3]
!P.CHARSIZE = 1.4
!X.TICKFORMAT = '(I6)'
!Y.TICKFORMAT = '(I5)'

sj[1] = '/del_quest.dat'
del = strjoin(sj)
openr, 1, del
S1 = fstat(1) 
rec_line = S1.SIZE/63.0 ;
print, rec_line
cal = dblarr(4, rec_line);  dblarr(col_no,row_no) 
readf, 1, cal
close, 1 

plot, cal(0,*), cal(1,*), /xstyle, /ystyle, $
		  title='TRUE-ESTI on LRS x-axis', $
		  background = 255, color = 0, $
		  symsize = 0.1, psym = 2, $
                  xrange=[ts, tf], $
;                  yrange=[-10, 10], $
                  ytitle='arcsec'

; xyouts, 0.90, 0.93, asn, /normal, align = 1, color = 0, charsize = 1.0 
; xyouts, 0.90, 0.93, year, /normal, align = 1, color = 0, charsize = 1.0
; xyouts, 0.94, 0.93, yday, /normal, align = 1, color = 0, charsize = 1.0

;for i=1, num_interval  do begin
;  xint = interval/3600.0 * i 
;  plots, [xint,xint], [-50, 50], linestyle = 1, color = 2
;endfor

plot, cal(0,*), cal(2,*),  /xstyle, /ystyle, $
		  title='TRUE-ESTI on LRS y-axis', $
		  background = 255, color = 0, $
		  symsize = 0.1, psym = 2, $
                  xrange=[ts, tf], $
;                  yrange=[-10, 10], $
                  ytitle='arcsec'

;for i=1, num_interval  do begin
;  xint = interval/3600.0 * i 
;  plots, [xint,xint], [-50, 50], linestyle = 1, color = 2
;endfor

plot, cal(0,*), cal(3,*), /xstyle, /ystyle, $
		  title='TRUE-ESTI on LRS z-axis', $
		  background = 255, color = 0, $
		  symsize = 0.1, psym = 2, $
                  xrange=[ts, tf], $
;                  yrange=[-100, 100], $
                  xtitle='time(hour)', $
                  ytitle='arcsec'

;for i=1, num_interval  do begin
;  xint = interval/3600.0 * i 
;  plots, [xint,xint], [-500, 500], linestyle = 1, color = 2
;endfor

sj[1] = '/quat_true_quest.png'
b1q = strjoin(sj)

image = tvrd()
loadct, 0
write_png, b1q, image

sj[1] = '/quat_true_quest.sta'
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

