import straxen
import saltax
import cutax
import strax
from immutabledict import immutabledict
from cutax.cut_lists.basic import BasicCuts
import warnings

# straxen XENONnT options/configuration
XNT_COMMON_OPTS = straxen.contexts.xnt_common_opts
XNT_COMMON_CONFIG = straxen.contexts.xnt_common_config
XNT_SIMULATION_CONFIG = straxen.contexts.xnt_simulation_config

# saltax options overrides
SXNT_COMMON_OPTS_REGISTER = XNT_COMMON_OPTS['register'].copy()
SXNT_COMMON_OPTS_REGISTER.remove(straxen.PulseProcessing)
SXNT_COMMON_OPTS_REGISTER = [saltax.SPulseProcessing] + SXNT_COMMON_OPTS_REGISTER
XNT_COMMON_OPTS_OVERRIDE = dict(
    register=SXNT_COMMON_OPTS_REGISTER,
)
SXNT_COMMON_OPTS = XNT_COMMON_OPTS.update(XNT_COMMON_OPTS_OVERRIDE)

# saltax configuration overrides
SCHANNEL_STARTS_AT = -494
XNT_COMMON_CONFIG_OVERRIDE = dict(
    channel_map=immutabledict(
        # (Minimum channel, maximum channel)
        # Channels must be listed in a ascending order!
        stpc=(SCHANNEL_STARTS_AT, SCHANNEL_STARTS_AT+493), # Salted TPC channels
        tpc=(0, 493),
        he=(500, 752),  # high energy
        aqmon=(790, 807),
        aqmon_nv=(808, 815),  # nveto acquisition monitor
        tpc_blank=(999, 999),
        mv=(1000, 1083),
        aux_mv=(1084, 1087),  # Aux mv channel 2 empty  1 pulser  and 1 GPS
        mv_blank=(1999, 1999),
        nveto=(2000, 2119),
        nveto_blank=(2999, 2999)
    ),
)
SXNT_COMMON_CONFIG = XNT_COMMON_CONFIG.update(XNT_COMMON_CONFIG_OVERRIDE)

# saltax modes supported
SALTAX_MODES = ['data', 'simu', 'salt']

# cutax XENONnT contexts
XENONNT_OFFLINE = cutax.contexts.xenonnt_offline
XENONNT_SIMULATION = cutax.contexts.xenonnt_sim_base
DEFAULT_XEDOCS_VERSION = cutax.contexts.DEFAULT_XEDOCS_VERSION


def xenonnt_salted(output_folder: str = './strax_data',
                   xedocs_version=DEFAULT_XEDOCS_VERSION,
                   cut_list=BasicCuts, auto_register=True,
                   faxconf_version="sr0_v4",
                   cmt_version="global_v11",
                   cmt_run_id="026000",
                   latest="sr0_v4",
                   **kwargs):
    # Based on cutax.xenonnt_sim_base()
    if faxconf_version != latest:
        warnings.warn(f'WARNING! You are using {faxconf_version.upper()} simulation context!'
                      f' Consider updating to a newer/pinned version ({latest.upper()})!')
    fax_conf='fax_config_nt_{:s}.json'.format(faxconf_version)

    # Based on straxen.contexts.xenonnt_online()
    context_options = dict(
        **SXNT_COMMON_OPTS,
        **kwargs,
    )
    context_config = dict(
        detector='XENONnT', # from straxen.contexts.xenonnt_simulation()
        fax_config=fax_conf, # from straxen.contexts.xenonnt_simulation()
        check_raw_record_overlaps=True, # from straxen.contexts.xenonnt_simulation()
        **SXNT_COMMON_CONFIG,
    )
    st = strax.Context(
        storage=strax.DataDirectory(output_folder)
        config=context_config,
        **context_options)
    st.register([straxen.DAQReader, saltax.SRawRecordsFromFaxNT])
    st.deregister_plugins_with_missing_dependencies()
        
    # Based on straxen.contexts.xenonnt()
    st.apply_cmt_version(cmt_version)
    if xedocs_version is not None:
        st.apply_xedocs_configs(version=xedocs_version, **kwargs)
    
    # Based on cutax.xenonnt_offline()
    # extra plugins to register
    st.set_config({'event_info_function': 'blinding_v11',
                   'avg_se_gain': 'bodega://se_gain?bodega_version=v1',
                   'g1': 'bodega://g1?bodega_version=v5',
                   'g2': 'bodega://g2?bodega_version=v5'})
    if auto_register:
        st.register_cuts()
    if cut_list is not None:
        st.register_cut_list(cut_list)

    # Based on straxen.xenonnt_simulation()
    _config_overlap=immutabledict(
            drift_time_gate='electron_drift_time_gate',
            drift_velocity_liquid='electron_drift_velocity',
            electron_lifetime_liquid='elife',
    )
    if straxen.utilix_is_configured(
            warning_message='Bad context as we cannot set CMT since we '
                            'have no database access'''):
        st.apply_cmt_version(cmt_version) 
    # Replace default cmt options with cmt_run_id tag + cmt run id
    cmt_options_full = straxen.get_corrections.get_cmt_options(st)
    # prune to just get the strax options
    cmt_options = {key: val['strax_option']
                   for key, val in cmt_options_full.items()}
    # First, fix gain model for simulation
    st.set_config({'gain_model_mc': 
                        ('cmt_run_id', cmt_run_id, *cmt_options['gain_model'])})
    fax_config_override_from_cmt = dict()
    for fax_field, cmt_field in _config_overlap.items():
        value = cmt_options[cmt_field]
        # URL configs need to be converted to the expected format
        if isinstance(value, str):
            opt_cfg = cmt_options_full[cmt_field]
            version = straxen.URLConfig.kwarg_from_url(value, 'version')
            # We now allow the cmt name to be different from the config name
            # WFSim expects the cmt name
            value = (opt_cfg['correction'], version, True)
        fax_config_override_from_cmt[fax_field] = ('cmt_run_id', cmt_run_id,
                                                   *value)
    st.set_config({'fax_config_override_from_cmt': fax_config_override_from_cmt})
    # and all other parameters for processing
    for option in cmt_options:
        value = cmt_options[option]
        if isinstance(value, str):
            # for URL configs we can just replace the run_id keyword argument
            # This will become the proper way to override the run_id for cmt configs
            st.config[option] = straxen.URLConfig.format_url_kwargs(value, run_id=cmt_run_id)
        else:
            # FIXME: Remove once all cmt configs are URLConfigs
            st.config[option] = ('cmt_run_id', cmt_run_id, *value)

    return st

    


def sxenonnt(saltax_mode,
             xedocs_version=DEFAULT_XEDOCS_VERSION,
             cut_list=BasicCuts, auto_register=True,
             faxconf_version="sr0_v4",
             cmt_version="global_v11",
             wfsim_registry='RawRecordsFromFaxNT',
             cmt_run_id="026000",
             latest="sr0_v4",
             **kwargs):
    assert saltax_mode in SALTAX_MODES, "saltax_mode must be one of %s"%(SALTAX_MODES)
    
    if saltax_mode == 'data':
        return XENONNT_OFFLINE(
            xedocs_version=xedocs_version,
            cut_list=cut_list,
            auto_register=auto_register,
            **kwargs)
    elif saltax_mode == 'simu':
        return XENONNT_SIMULATION(
            faxconf_version=faxconf_version,
            cmt_version=cmt_version,
            wfsim_registry=wfsim_registry,
            cmt_run_id=cmt_run_id,
            latest=latest,
            cut_list=cut_list,
            **kwargs)
    