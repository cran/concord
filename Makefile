#
# ${R_HOME}/src/library/concord/Makefile

srcdir = .
top_srcdir = ../../..
top_builddir = ../../..
subdir = src/library/concord

include $(top_builddir)/Makeconf

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = INDEX TITLE DESCRIPTION.in Makefile.in

pkg = concord

all: Makefile DESCRIPTION
	@echo "building package '$(pkg)'"
	@$(MKINSTALLDIRS) $(top_builddir)/library/$(pkg)/R
	@(f=$${TMPDIR:-/tmp}/R$$$$; \
	  cat `ls $(srcdir)/R/*.R` > $${f}; \
	  $(top_srcdir)/tools/move-if-change $${f} \
	    $(top_builddir)/library/$(pkg)/R/$(pkg))
	@for f in INDEX TITLE; do \
	if test -f $(srcdir)/$${f}; then \
	  $(INSTALL_DATA) $(srcdir)/$${f} \
	    $(top_builddir)/library/$(pkg); \
	fi; \
	done
	@if test -f $(top_builddir)/$(subdir)/DESCRIPTION; then \
	  $(INSTALL_DATA) $(top_builddir)/$(subdir)/DESCRIPTION \
	    $(top_builddir)/library/$(pkg); \
	  echo "Built: R" @VERSION@\;  @R_PLATFORM@\;  `date` >> \
	    $(top_builddir)/library/$(pkg)/DESCRIPTION; \
	fi
	@if test -d $(srcdir)/data; then \
	  $(MKINSTALLDIRS) $(top_builddir)/library/$(pkg)/data; \
	  for f in `ls -d $(srcdir)/data/* | sed '/CVS/d'`; do \
	    $(INSTALL_DATA) $${f} $(top_builddir)/library/$(pkg)/data; \
	  done; \
	fi
	@if test -d $(srcdir)/man; then \
	  $(MKINSTALLDIRS) $(top_builddir)/library/$(pkg)/man; \
	  (f=$${TMPDIR:-/tmp}/R$$$$; \
	    cat $(srcdir)/man/*.Rd > $${f}; \
	    $(top_srcdir)/tools/move-if-change $${f} \
	      $(top_builddir)/library/$(pkg)/man/$(pkg).Rd); \
	fi
	@if test -d src; then \
	  (cd src && $(MAKE)) || exit 1; \
	fi

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
DESCRIPTION: $(srcdir)/DESCRIPTION.in $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

mostlyclean: clean
clean:
	@if test -d src; then (cd src && $(MAKE) $@); fi
distclean:
	@if test -d src; then (cd src && $(MAKE) $@); fi
	-@rm -f Makefile DESCRIPTION
maintainer-clean: distclean

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done
	@for d in R data exec man src; do \
	  if test -d $(srcdir)/$${d}; then \
	    ((cd $(srcdir); \
	          $(TAR) -c -f - $(DISTDIR_TAR_EXCLUDE) $${d}) \
	        | (cd $(distdir); $(TAR) -x -f -)) \
	      || exit 1; \
	  fi; \
	done