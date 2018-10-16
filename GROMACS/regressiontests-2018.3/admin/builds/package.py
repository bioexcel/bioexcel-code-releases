def do_build(context):
    version = context.params.get('PACKAGE_VERSION_STRING', Parameter.string)
    release = context.params.get('RELEASE', Parameter.bool)
    if not version:
        version = 'unknown'
    if not release:
        version += '-dev'

    package_name = 'regressiontests-' + version

    context.make_archive(package_name, use_git=True, prefix=package_name)

    package_name += '.tar.gz'
    context.write_package_info(Project.REGRESSIONTESTS, package_name, version)
