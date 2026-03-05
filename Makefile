SHELL := /bin/bash

# SET THESE
REMOTE_USER ?=
REMOTE_HOST ?=koa.its.hawaii.edu



REMOTE := $(REMOTE_USER)@$(REMOTE_HOST)
SSH_ALIAS ?= koa
SSH_CONFIG ?= $(HOME)/.ssh/config
SSH_PERSIST ?= 8h
SSH_TARGET ?= $(SSH_ALIAS)

# Default output directories when OUTDIR is not provided.
REMOTE_HOME ?= ~
LOCAL_HOME ?= $(HOME)
GOAL_ITEMS := $(filter-out upload download,$(MAKECMDGOALS))
EFFECTIVE_ITEMS := $(strip $(if $(ITEMS),$(ITEMS),$(GOAL_ITEMS)))

.PHONY: help login upload download setup-ssh-cache cache-open cache-close cache-check

ifneq (,$(filter upload download,$(MAKECMDGOALS)))
$(GOAL_ITEMS):
	@:
endif

help:
	@echo "High performance computing helper functions (all functions run from local computer):"
	@echo "  make login   (quick login)"
	@echo ""
	@echo "  make cache-open   (open/reuse master ssh session)"
	@echo "  make cache-check  (check whether master session exists)"
	@echo "  make cache-close  (close master ssh session)"
	@echo "  make setup-ssh-cache [SSH_ALIAS=koa] [REMOTE_USER=...] [REMOTE_HOST=...] [SSH_PERSIST=8h]"
	@echo ""
	@echo "  make upload file1 dir1 [OUTDIR=~/remote/path/]"
	@echo "  make download file1 dir1 [OUTDIR=~/local/path/]"
	@echo ""
	@echo "Set REMOTE_USER and REMOTE_HOST in the make file or pass at runtime or use setup-ssh-cache to store them."

login:
	ssh $(SSH_TARGET)

upload:
	@[ -n "$(EFFECTIVE_ITEMS)" ] || (echo "Usage: make upload ITEMS=\"file1 dir1\" [OUTDIR=<path>]"; echo "   or: make upload file1 dir1 [OUTDIR=<path>]"; exit 1)
	@if ! ssh -O check "$(SSH_TARGET)" >/dev/null 2>&1; then \
		echo "Opening SSH cache for $(SSH_TARGET) (Duo may be required once)..."; \
		ssh -MNf "$(SSH_TARGET)"; \
	fi
	@dest="$(if $(OUTDIR),$(OUTDIR),$(REMOTE_HOME))"; \
	scp -r $(EFFECTIVE_ITEMS) "$(SSH_TARGET):$$dest"

download:
	@[ -n "$(EFFECTIVE_ITEMS)" ] || (echo "Usage: make download ITEMS=\"~/remote/file ~/remote/dir\" [OUTDIR=<path>]"; echo "   or: make download ~/remote/file ~/remote/dir [OUTDIR=<path>]"; exit 1)
	@if ! ssh -O check "$(SSH_TARGET)" >/dev/null 2>&1; then \
		echo "Opening SSH cache for $(SSH_TARGET) (Duo may be required once)..."; \
		ssh -MNf "$(SSH_TARGET)"; \
	fi
	@dest="$(if $(OUTDIR),$(OUTDIR),$(LOCAL_HOME))"; \
	set --; \
	for item in $(EFFECTIVE_ITEMS); do \
		set -- "$$@" "$(SSH_TARGET):$$item"; \
	done; \
	scp -r "$$@" "$$dest"

cache-open:
	@if ssh -O check "$(SSH_TARGET)" >/dev/null 2>&1; then \
		echo "SSH cache already open for $(SSH_TARGET)"; \
	else \
		echo "Opening SSH cache for $(SSH_TARGET)..."; \
		ssh -MNf "$(SSH_TARGET)"; \
	fi

cache-check:
	@if ssh -O check "$(SSH_TARGET)" >/dev/null 2>&1; then \
		echo "SSH cache is open for $(SSH_TARGET)"; \
	else \
		echo "SSH cache is not open for $(SSH_TARGET)"; \
		exit 1; \
	fi

cache-close:
	@if ssh -O check "$(SSH_TARGET)" >/dev/null 2>&1; then \
		ssh -O exit "$(SSH_TARGET)"; \
		echo "Closed SSH cache for $(SSH_TARGET)"; \
	else \
		echo "No SSH cache session to close for $(SSH_TARGET)"; \
	fi

setup-ssh-cache:
	@[ -n "$(REMOTE_USER)" ] || (echo "Missing remote user. Set REMOTE_USER in Makefile or pass REMOTE_USER=<user>"; exit 1)
	@[ -n "$(REMOTE_HOST)" ] || (echo "Missing remote host. Set REMOTE_HOST in Makefile or pass REMOTE_HOST=<host>"; exit 1)
	@mkdir -p "$$(dirname "$(SSH_CONFIG)")"
	@touch "$(SSH_CONFIG)"
	@tmp="$$(mktemp)"; \
	start="# >>> ssh $(SSH_ALIAS) >>>"; \
	end="# <<< ssh $(SSH_ALIAS) <<<"; \
	awk -v s="$$start" -v e="$$end" '\
		$$0==s {inblk=1; next} \
		$$0==e {inblk=0; next} \
		!inblk {print}\
	' "$(SSH_CONFIG)" > "$$tmp"; \
	printf "%s\n" "$$start" >> "$$tmp"; \
	printf "Host %s %s\n" "$(SSH_ALIAS)" "$(REMOTE_HOST)" >> "$$tmp"; \
	printf "  HostName %s\n" "$(REMOTE_HOST)" >> "$$tmp"; \
	printf "  User %s\n" "$(REMOTE_USER)" >> "$$tmp"; \
	printf "  ControlMaster auto\n" >> "$$tmp"; \
	printf "  ControlPath ~/.ssh/cm-%%r@%%h:%%p\n" >> "$$tmp"; \
	printf "  ControlPersist %s\n" "$(SSH_PERSIST)" >> "$$tmp"; \
	printf "%s\n" "$$end" >> "$$tmp"; \
	mv "$$tmp" "$(SSH_CONFIG)"; \
	chmod 600 "$(SSH_CONFIG)"; \
	echo "Updated $(SSH_CONFIG) with Host $(SSH_ALIAS)"
