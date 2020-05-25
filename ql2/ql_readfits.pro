function ql_readfits, filename, hd, EXTEN_NO=extension

if (arg_present(extension) ne 1) then extension=0

; set up a catch to make sure there isn't a problem when
; trying to open a file
catch, error_status

; This statement begins the error handler:
if error_status ne 0 then begin
    print, 'Error index: ', error_status
    print, 'Error message: ', !ERROR_STATE.MSG
endif

; Check to see if the file is completely written, and proceed accordingly
if (error_status eq 0) then begin
   fits_read = readfits(filename, hd, EXTEN_NO = extension)
   if n_elements(fits_read) eq 1 then begin
      ;; if there is no data in extension 0, try extension 1
      ;; (e.g. NIFS data)
      print, 'No data in extension ', extension, ' reading extension 1'
      fits_read =  readfits(filename, hd1, EXTEN_NO = 1)
   endif
   
    return, fits_read
endif else begin
    print, 'There was an error reading the fits file ', filename, '.'
    return, -1
endelse

end
